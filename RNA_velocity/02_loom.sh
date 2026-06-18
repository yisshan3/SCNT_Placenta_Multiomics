#!/bin/sh

# ==============================================================================
# Description: Generate loom files containing spliced and unspliced reads from Cell Ranger Arc output using velocyto.
# ==============================================================================


# Setup Environment ------------------------------------------------------------
cd /home4/ssyi/Mouse_Placenta/codes/4.RNA_velocity || exit
mkdir -p loom

# Activate conda environment
conda activate scVelo


# Prepare Input Metadata -------------------------------------------------------
touch input_metadata.txt
for sample in WT_batch1_125_PWK WT_batch2_145_PWK SCNT_batch3_125_PWK WT_batch4_185_PWK SCNT_batch5_145_PWK SCNT_batch6_185_PWK; do
    echo -e "/home4/ssyi/Mouse_Placenta/2.cellranger-arc/${sample}\t${sample%_PWK}" >> input_metadata.txt
done


# Generate loom ----------------------------------------------------------------
# Define reference files
rmsk_gtf="/home1/ssyi/annotate/mm/Repeats/mm10_repeatMasker_allTracks.gtf"
cellranger_gtf="/home1/ssyi/annotate/mm/SNP/cellranger-arc/mm10_PWK_PhJ_N-masked/genes/genes.gtf"

cat input_metadata.txt|while read i;do
    arr=($i)   
    dir=${arr[0]} 
    sample=${arr[1]}
    
    nohup velocyto run -m $rmsk_gtf ${dir}/outs/gex_possorted_bam.bam $cellranger_gtf \
        -@ 10 --samtools-memory 5000 \
        -e ${sample} \
        -b ${dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
        -o /home4/ssyi/Mouse_Placenta/codes/4.RNA_velocity/loom/ \
        > ./${sample}_loom.log 2>&1 &
done

