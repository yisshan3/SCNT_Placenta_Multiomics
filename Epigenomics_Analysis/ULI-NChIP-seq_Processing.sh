#!/bin/bash

# ==============================================================================
# Description: Upstream processing pipeline for H3K27me3 ChIP-seq data, inluding QC, trimming, alignment, deduplication, bigwig generation.
# ==============================================================================


# Setup ------------------------------------------------------------------------
E75EXE_H3K27me3="/home4/ssyi/Mouse_Placenta/Public_data/WMZ_E75EXE_H3K27me3NChIP"
bowtie2_index="/home/share/bowtie2_index/mm10"

cd "${E75EXE_H3K27me3}" || exit


# FastQC & MultiQC -------------------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/1.qc"
cd "${E75EXE_H3K27me3}/1.qc" || exit

for id in "${E75EXE_H3K27me3}/0.raw/"*.fq.gz; do
    fastqc -t 10 "${id}" -o "${E75EXE_H3K27me3}/1.qc" &
done
wait 

multiqc -f "${E75EXE_H3K27me3}/1.qc" -o "${E75EXE_H3K27me3}/1.qc" -n raw_multiqc_report


# Trim Galore ------------------------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/2.trim_galore"
cd "${E75EXE_H3K27me3}/2.trim_galore" || exit

for i in "${E75EXE_H3K27me3}/0.raw/"*R1.fq.gz; do
    r1="${i}"
    r2="${i/R1.fq.gz/R2.fq.gz}"
    trim_galore -j 10 -q 25 --length 50 --trim-n --paired --fastqc \
        -o "${E75EXE_H3K27me3}/2.trim_galore" "${r1}" "${r2}" &
done
wait

multiqc "${E75EXE_H3K27me3}/2.trim_galore" -o "${E75EXE_H3K27me3}/2.trim_galore" -n clean_multiqc_report


# Bowtie2 Alignment ------------------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/3.bowtie2"
cd "${E75EXE_H3K27me3}/3.bowtie2" || exit

for i in "${E75EXE_H3K27me3}/2.trim_galore/"*.R1_val_1.fq.gz; do
    sample=$(basename "${i}")
    sample="${sample%%.*}"
    r1="${i}"
    r2="${i/R1_val_1.fq.gz/R2_val_2.fq.gz}"
    
    bowtie2 -p 10 --end-to-end --very-sensitive \
        --no-mixed --no-discordant --no-unal -x "${bowtie2_index}" \
        -1 "${r1}" -2 "${r2}" -S "${E75EXE_H3K27me3}/3.bowtie2/${sample}.sam" \
        > "${E75EXE_H3K27me3}/3.bowtie2/${sample}.log" 2>&1 &
done
wait

# Parse Bowtie2 alignment logs
for id in *.log; do
    sample=$(basename "${id}" ".log")
    clean_pairs=$(grep "paired" "${id}" | sed 's/^ *//' | cut -d " " -f 1)
    aligned_rate=$(grep "overall" "${id}" | cut -d " " -f 1 | tr "\n" "\t" | sed "s/\t$/\n/")
    one=$(grep "concordantly" "${id}" | grep "exactly" | sed 's/^ *//' | cut -d " " -f 1)
    one_rate=$(grep "concordantly" "${id}" | grep "exactly" | sed 's/^ *//' | cut -d "(" -f 2 | cut -d ")" -f 1)
    onemore=$(grep "concordantly" "${id}" | grep ">" | sed 's/^ *//' | cut -d " " -f 1)
    onemore_rate=$(grep "concordantly" "${id}" | grep ">" | sed 's/^ *//' | cut -d "(" -f 2 | cut -d ")" -f 1)
    aligned_pairs=$(expr "${one}" + "${onemore}")
    echo -e "${sample}\t${clean_pairs}\t${one}\t${one_rate}\t${onemore}\t${onemore_rate}\t${aligned_pairs}\t${aligned_rate}" >> bt2_align_stat.txt
done
sed -i '1s/^/sample\ttotal_pairs\tproperly_paired\tratio\tmore_1_time\tratio\taligned_pairs\taligned_rate\n/' bt2_align_stat.txt


# Samtools Sort & Index --------------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/4.samtools"
cd "${E75EXE_H3K27me3}/4.samtools" || exit

for i in "${E75EXE_H3K27me3}/3.bowtie2/"*.sam; do
    sample=$(basename "${i}" ".sam")
    samtools sort -@ 8 -O bam -o "${E75EXE_H3K27me3}/4.samtools/${sample}.sorted.bam" "${i}" &
done
wait

for i in *.sorted.bam; do
    samtools index -@ 10 "${i}" &
done
wait


# Sambamba Markdup & Merge -----------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/5.sambamba/merge"
cd "${E75EXE_H3K27me3}/5.sambamba" || exit

# Remove duplicates
for i in "${E75EXE_H3K27me3}/4.samtools/"*.sorted.bam; do
    sample=$(basename "${i}" ".sorted.bam")
    sambamba markdup -r -t 8 -p "${i}" "${E75EXE_H3K27me3}/5.sambamba/${sample}.rmdup.bam" \
        > "${E75EXE_H3K27me3}/5.sambamba/${sample}.rmdup.log" 2>&1 &
done
wait

# Flagstat
for i in "${E75EXE_H3K27me3}/5.sambamba/"*.rmdup.bam; do
    sample=$(basename "${i}" ".rmdup.bam")
    samtools flagstat -@ 5 "${i}" > "${sample}.flagstat" &
done
wait

for i in *.flagstat; do
    sample=$(basename "${i}" ".flagstat")
    dedup_pairs=$(grep "read1" "${i}" | cut -d ' ' -f1)
    echo -e "${sample}\t${dedup_pairs}" >> flagstat.txt
done
sed -i '1s/^/sample\tdedup_pairs\n/' flagstat.txt

ln -s ../3.bowtie2/bt2_align_stat.txt ./
awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a){print $0,a[$1]}}' flagstat.txt bt2_align_stat.txt > bt2_align_stat_sambamba.txt

# Merge replicates
for k in $(for i in "${E75EXE_H3K27me3}/5.sambamba/"*.rmdup.bam; do o=$(basename "$i"); echo "${o%%_rep*}"; done | sort | uniq); do 
    samtools merge -@ 8 "${E75EXE_H3K27me3}/5.sambamba/merge/${k}.rmdup.bam" "${E75EXE_H3K27me3}/5.sambamba/${k}"*rmdup.bam & 
done
wait

for i in "${E75EXE_H3K27me3}/5.sambamba/merge/"*.rmdup.bam; do
    samtools index -@ 10 "${i}" &
done
wait


# DeepTools --------------------------------------------------------------------
mkdir -p "${E75EXE_H3K27me3}/6.deeptools/bamCoverage/rep"
mkdir -p "${E75EXE_H3K27me3}/6.deeptools/bamCoverage/merge"

# Replicates coverage
for i in "${E75EXE_H3K27me3}/5.sambamba/"*.rmdup.bam; do
    sample=$(basename "${i}" ".rmdup.bam")
    bamCoverage -p 10 -b "${i}" -bs 50 --normalizeUsing RPKM \
        -of "bigwig" -o "${E75EXE_H3K27me3}/6.deeptools/bamCoverage/rep/${sample}.bw" &
done

# Merged coverage
for i in "${E75EXE_H3K27me3}/5.sambamba/merge/"*.rmdup.bam; do
    sample=$(basename "${i}" ".rmdup.bam")
    bamCoverage -p 10 -b "${i}" -bs 50 --normalizeUsing RPKM \
        -of "bigwig" -o "${E75EXE_H3K27me3}/6.deeptools/bamCoverage/merge/${sample}.bw" &
done
wait

