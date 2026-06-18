#!/bin/bash

# ==============================================================================
# Description: Pipeline for single-nuclei SNPsplit. Extracts barcodes, splits BAM files using sinto, converts to FASTQ, aligns RNA/ATAC reads to an N-masked genome, runs SNPsplit, and calculates allele ratios.
# ==============================================================================

set -euo pipefail

# Utility Functions-------------------------------------------------------------
function help() { 
    cat << EOF
Usage: sh 03_SingleNuclei_SNPsplit.sh -s <sample>

Required:
  -s  Sample name, including:
      WT_batch1_125_PWK    WT_batch2_145_PWK    SCNT_batch3_125_PWK
      WT_batch4_185_PWK    SCNT_batch5_145_PWK  SCNT_batch6_185_PWK

Options:
  -h, --help    Show this help message
EOF
}

function exit1(){ 
    echo "ERROR: $1" >&2
    exit 1
}


# Argument Parsing--------------------------------------------------------------
# Transform long options to short ones
for arg in "$@"; do
    shift
    case "$arg" in
        "--help")   set -- "$@" "-h" ;;
        *)        set -- "$@" "$arg"
    esac
done

# Parse command-line arguments
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;;
    -s) sample=$2; shift 2;;
     *) echo "error: no such option $1. -h|--help for help";exit 1;;
esac
done

# Validate required arguments
if [ -z "${sample:-}" ]; then
    exit1 "Sample name is required. Use -h for help."
fi


# Directory Setup --------------------------------------------------------------
ref_dir="/home1/ssyi/annotate/mm/SNP/PWK_PhJ_Single_strain"
base_dir="/home4/ssyi/Mouse_Placenta/2.cellranger-arc/${sample}/outs"
snp_dir="${base_dir}/SNP"

barcode_dir="${snp_dir}/1.whitelist_barcode"
sinto_dir="${snp_dir}/2.sinto_split"
fq_dir="${snp_dir}/3.bam2fastq"
align_dir="${snp_dir}/4.alignment"
snpsplit_dir="${snp_dir}/5.SNPsplit"
ratio_dir="${snp_dir}/6.calculate_allele_ratio"


# Run pipeline------------------------------------------------------------------
echo "========================================================================"
echo "Starting SNP allele-specific pipeline for sample: ${sample}"
echo "========================================================================"

# --------------------------------------------------------------------------
# 01. Filter Whitelist Barcodes
# --------------------------------------------------------------------------
echo "--- 01. Filtering Whitelist Barcodes ---"
mkdir -p "${barcode_dir}"
cd "${base_dir}" || exit1 "Cannot access ${base_dir}"

if [ ! -f "./filtered_feature_bc_matrix/barcodes.tsv.gz" ]; then
    exit1 "Barcode file not found: ./filtered_feature_bc_matrix/barcodes.tsv.gz"
fi

# Extract GEX barcodes and remove the '-1' suffix
zcat ./filtered_feature_bc_matrix/barcodes.tsv.gz | cut -d "-" -f 1 > "${barcode_dir}/whitelist_barcode.tab"


# --------------------------------------------------------------------------
# 02. Sinto Split BAM by Barcode
# --------------------------------------------------------------------------
echo "--- 02. Splitting BAM files using Sinto ---"
mkdir -p "${sinto_dir}/RNA" "${sinto_dir}/ATAC"
cd "${barcode_dir}" || exit1 "Cannot access ${barcode_dir}"

# Reformat barcodes for Sinto
awk 'BEGIN{FS=OFS="\t"}{print $1"-1", $1"-1"}' whitelist_barcode.tab > barcodes.tab

# Verify input BAM files
if [ ! -f "${base_dir}/gex_possorted_bam.bam" ]; then
    exit1 "RNA BAM not found: ${base_dir}/gex_possorted_bam.bam"
fi
if [ ! -f "${base_dir}/atac_possorted_bam.bam" ]; then
    exit1 "ATAC BAM not found: ${base_dir}/atac_possorted_bam.bam"
fi

# Split RNA and ATAC BAMs
sinto filterbarcodes -b "${base_dir}/gex_possorted_bam.bam" -c barcodes.tab -p 8 --barcodetag "CB" --outdir "${sinto_dir}/RNA" > "${sinto_dir}/RNA/sinto_rna.log" 2>&1 &
sinto filterbarcodes -b "${base_dir}/atac_possorted_bam.bam" -c barcodes.tab -p 8 --barcodetag "CB" --outdir "${sinto_dir}/ATAC" > "${sinto_dir}/ATAC/sinto_atac.log" 2>&1 &
wait

# Index BAM files
for i in "${sinto_dir}/RNA/"*.bam "${sinto_dir}/ATAC/"*.bam; do
    samtools index -@ 10 "${i}" &
done
wait


# --------------------------------------------------------------------------
# 03. BAM to FASTQ Conversion
# --------------------------------------------------------------------------
echo "--- 03. Converting BAM to FASTQ ---"
mkdir -p "${fq_dir}/ATAC" "${fq_dir}/RNA"

# ATAC
cd "${fq_dir}/ATAC" || exit1 "Cannot access ${fq_dir}/ATAC"
for i in "${sinto_dir}/ATAC/"*.bam; do
    sample_name=$(basename "${i}" ".bam")
    samtools sort -n -@ 8 -o "./${sample_name}.nsort.bam" "${i}" &
done
wait

for i in *.nsort.bam; do
    sample_name=$(basename "${i}" ".nsort.bam")
    samtools fastq -1 "${sample_name}.R1.fq.gz" -2 "${sample_name}.R2.fq.gz" -0 /dev/null -s /dev/null -n -@ 8 "${i}" &
done
wait
rm -f *.nsort.bam

# RNA
cd "${fq_dir}/RNA" || exit1 "Cannot access ${fq_dir}/RNA"
for i in "${sinto_dir}/RNA/"*.bam; do
    sample_name=$(basename "${i}" ".bam")
    bash -c "samtools bam2fq -@ 10 ${i} 2> ${sample_name}.log | gzip > ${sample_name}.fq.gz" & 
done
wait


# --------------------------------------------------------------------------
# 04. Alignment to N-masked Genome
# --------------------------------------------------------------------------
echo "--- 04. Alignment ---"
mkdir -p "${align_dir}/ATAC" "${align_dir}/RNA"

# ATAC (Bowtie2)
cd "${align_dir}/ATAC" || exit1 "Cannot access ${align_dir}/ATAC"
for i in "${fq_dir}/ATAC/"*.R1.fq.gz; do
    r1="${i}"
    r2="${i/R1.fq/R2.fq}"
    sample_name=$(basename "${i}" ".R1.fq.gz")

    bowtie2 -p 10 -x "${ref_dir}/PWK_PhJ_N-masked_bowtie2_index/PWK_PhJ_N-masked" \
        --end-to-end --very-sensitive -X 2000 \
        --no-mixed --no-discordant --no-unal \
        --time -1 "${r1}" -2 "${r2}" -S "${sample_name}.sam" > "${sample_name}.bt2.log" 2>&1 &
done
wait

# Parse Bowtie2 alignment rates
for id in *.bt2.log; do
    sample_name=$(basename "${id}" ".bt2.log")
    rate=$(grep "overall" "${id}" | cut -d " " -f 1 | tr "\n" "\t" | sed "s/\t$/\n/")
    echo -e "${sample_name}\t${rate}" >> bt2_align_rate.txt
done
sed -i '1s/^/sample\taligned_rate\n/' bt2_align_rate.txt

for i in *.sam; do
    sample_name=$(basename "${i}" ".sam")
    samtools view -@ 10 -bS "${i}" > "./${sample_name}.bam" &
done
wait
rm -f *.sam

# RNA (HISAT2)
cd "${align_dir}/RNA" || exit1 "Cannot access ${align_dir}/RNA"
for i in "${fq_dir}/RNA/"*.fq.gz; do
    sample_name=$(basename "${i}" ".fq.gz")
    hisat2 -p 10 -x "${ref_dir}/PWK_PhJ_N-masked_hisat2_index/PWK_PhJ_N-masked" \
        --dta-cufflinks --no-unal --no-softclip -k 1 \
        -t -U "${i}" -S "${sample_name}.sam" > "${sample_name}.ht2.log" 2>&1 &
done
wait

# Parse HISAT2 alignment rates
for id in *.ht2.log; do
    sample_name=$(basename "${id}" ".ht2.log")
    rate=$(grep "overall" "${id}" | cut -d " " -f 1 | tr "\n" "\t" | sed "s/\t$/\n/")
    echo -e "${sample_name}\t${rate}" >> ht2_align_rate.txt
done
sed -i '1s/^/sample\taligned_rate\n/' ht2_align_rate.txt

for i in *.sam; do
    sample_name=$(basename "${i}" ".sam")
    samtools view -@ 10 -bS "${i}" > "./${sample_name}.bam" &
done
wait
rm -f *.sam


# --------------------------------------------------------------------------
# 05. SNPsplit
# --------------------------------------------------------------------------
echo "--- 05. SNPsplit ---"
mkdir -p "${snpsplit_dir}/ATAC" "${snpsplit_dir}/RNA"
snp_list="${ref_dir}/all_SNPs_PWK_PhJ_GRCm38_modified.txt"
if [ ! -f "${snp_list}" ]; then
    exit1 "SNP list not found: ${snp_list}"
fi

# Split ATAC reads by allele
cd "${snpsplit_dir}/ATAC" || exit1 "Cannot access ${snpsplit_dir}/ATAC"
for i in "${align_dir}/ATAC/"*.bam; do
    sample_name=$(basename "${i}" ".bam")
    SNPsplit --paired --conflicting -o "./" --snp_file "${snp_list}" "${i}" > "./${sample_name}.log" 2>&1 &
done
wait

# Split RNA reads by allele
cd "${snpsplit_dir}/RNA" || exit1 "Cannot access ${snpsplit_dir}/RNA"
for i in "${align_dir}/RNA/"*.bam; do
    sample_name=$(basename "${i}" ".bam")
    SNPsplit --single_end --conflicting -o "./" --snp_file "${snp_list}" "${i}" > "./${sample_name}.log" 2>&1 &
done
wait


# --------------------------------------------------------------------------
# 06. Calculate Allele Ratio
# --------------------------------------------------------------------------
echo "--- 06. Calculating Allele Ratios ---"
mkdir -p "${ratio_dir}"

cd "${snpsplit_dir}/ATAC" || exit1 "Cannot access ${snpsplit_dir}/ATAC"
summary_file="${ratio_dir}/ATAC_SNPsplit_summary.tab"
> "${summary_file}"

for i in $(ls */*.SNPsplit_report.txt | grep -v ".sortedByName"); do
    CB=$(basename $i ".SNPsplit_report.txt")
    genome1_reads=$(grep "genome 1" $i | grep "%" | grep -v ":" | cut -d " " -f1)
    genome1_ratio=$(grep "genome 1" $i | grep "%" | grep -v ":" | cut -d "(" -f 2 | cut -d ")" -f 1)
    genome2_reads=$(grep "genome 2" $i | grep "%" | grep -v ":" | cut -d " " -f1)
    genome2_ratio=$(grep "genome 2" $i | grep "%" | grep -v ":" | cut -d "(" -f 2 | cut -d ")" -f 1)
    echo -e ${CB}"\t"${genome1_reads}"\t"${genome1_ratio}"\t"${genome2_reads}"\t"${genome2_ratio} >> "${summary_file}"
done
sed -i '1s/^/CB\tgenome1_reads\tgenome1_ratio\tgenome2_reads\tgenome2_ratio\n/' "${summary_file}"

echo "Finished processing sample: ${sample}"

echo "Summary file: ${ratio_dir}/ATAC_SNPsplit_summary.tab"


