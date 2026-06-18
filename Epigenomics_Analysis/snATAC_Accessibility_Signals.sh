#!/bin/bash

# ==============================================================================
# Description: Split ATAC-seq BAM files by cell barcodes using sinto, merge them by cell type, generate BigWig files, and compute accessibility signals over DARs.
# ==============================================================================

set -euo pipefail

# Setup Directories ------------------------------------------------------------
sinto_dir="/home4/ssyi/Mouse_Placenta/codes/9.sinto"
cd "${sinto_dir}" || exit


# Split BAM files by barcodes for each batch -----------------------------------
mkdir -p ./ATAC/Batch
Batch_dir="${sinto_dir}/ATAC/Batch"
cd "${Batch_dir}" || exit

samples=(
  "SCNT_batch3_125_PWK"
  "SCNT_batch5_145_PWK"
  "SCNT_batch6_185_PWK"
  "WT_batch1_125_PWK"
  "WT_batch2_145_PWK"
  "WT_batch4_185_PWK"
)

for i in "${samples[@]}"; do
  sample_dir="${Batch_dir}/${i}"
  mkdir -p "${sample_dir}/sinto"
  cd "${sample_dir}" || exit
  
  sinto filterbarcodes \
    -b "/home4/ssyi/Mouse_Placenta/2.cellranger-arc/${i}/outs/atac_possorted_bam.bam" \
    -c Cell_meta.txt -p 12 --barcodetag "CB" --outdir ./sinto > sinto.log 2>&1 &
done
wait


# Merge BAM files by Cell Type -------------------------------------------------
mkdir -p "${sinto_dir}/ATAC/Celltype"
cd "${sinto_dir}/ATAC/Celltype" || exit


# Merge WT
cat Tropho_celltypes.txt | while read -r i; do
    samtools merge --threads 6 -o "WT_${i}.bam" "${sinto_dir}/ATAC/Batch/"*WT*"/sinto/${i}.bam"
done
wait

# Merge SCNT
cat Tropho_celltypes.txt | while read -r i; do
    samtools merge --threads 6 -o "SCNT_${i}.bam" "${sinto_dir}/ATAC/Batch/"*SCNT*"/sinto/${i}.bam" &
done
wait


# Index and convert to BigWig --------------------------------------------------
# Index BAMs
for i in *.bam; do
    samtools index -@ 10 "${i}" &
done
wait

# Convert to bw
for i in ./*.bam; do
    prefix=$(basename "${i}" ".bam")
    bamCoverage -p 10 -b "${i}" -bs 50 --normalizeUsing RPKM -of "bigwig" -o "./${prefix}.bw" &
done
wait


# Compute Signals over DARs ----------------------------------------------------
DAR_dir="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Analysis/DARs"
cd "${DAR_dir}" || exit

# Compute Matrix
cat "${sinto_dir}/ATAC/Celltype/Precursor_celltypes.txt" | while read -r i; do
  computeMatrix scale-regions \
    -p 20 \
    -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
    -R "${DAR_dir}/SCNT_${i}_DARs.bed" \
    -S "${sinto_dir}/ATAC/Celltype/WT_${i}.bw" "${sinto_dir}/ATAC/Celltype/SCNT_${i}.bw" \
    --samplesLabel WT SCNT \
    --skipZeros \
    -o "SCNT_${i}_ATAC.gz" \
    --outFileSortedRegions "SCNT_${i}_ATAC.bed" &
done
wait

# Plot Heatmap
cat "${sinto_dir}/ATAC/Celltype/Precursor_celltypes.txt" | while read -r i; do
  plotHeatmap \
    --matrixFile "SCNT_${i}_ATAC.gz"  \
    --outFileName "SCNT_${i}_ATAC_heatmap.pdf" \
    --plotTitle "WT vs SCNT ${i}" \
    --whatToShow "plot, heatmap and colorbar" \
    --colorMap Blues \
    --perGroup \
    --heatmapHeight 30 \
    --heatmapWidth 3 &
done
wait

