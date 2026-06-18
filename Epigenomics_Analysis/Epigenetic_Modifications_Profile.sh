#!/bin/bash

# ==============================================================================
# Description: Calculate and plot deepTools profiles (Accessibility, H3K4me3, H3K27me3, H3K9me3) around DDARs and random regions.
# ==============================================================================


# Setup Directories ------------------------------------------------------------
cd /home4/ssyi/Mouse_Placenta/Figures || exit
DAR_dir="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Analysis/DARs"


# Accessibility ----------------------------------------------------------------
WT_Acc_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/WT/Accessibility"
SCNT_Acc_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/SCNT/Accessibility"   

# GV/CC DNase & GV/M2 ATAC
nohup computeMatrix scale-regions -p 20 -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
  -R "${DAR_dir}/DDARs.bed" "${DAR_dir}/Random_regions.bed" \
  -S "${WT_Acc_bw}/Dseq_Oocyte.bw" "${SCNT_Acc_bw}/cc_dnasei.bw" "${WT_Acc_bw}/GV_ATAC.bw" "${WT_Acc_bw}/MII_ATAC.bw" \
  --samplesLabel GV_DNase CC_DNase GV_ATAC MII_ATAC \
  --skipZeros \
  --outFileName GV_M2_CC_chrom_acc_signal_around_DDARs.gz \
  --outFileNameMatrix GV_M2_CC_chrom_acc_signal_around_DDARs.txt \
  --outFileSortedRegions GV_M2_CC_chrom_acc_signal_around_DDARs.bed &

nohup plotProfile -m GV_M2_CC_chrom_acc_signal_around_DDARs.gz \
  -out GV_M2_CC_chrom_acc_signal_around_DDARs.pdf \
  --colors '#a6cee3' '#1f78b4' '#b2df8a' '#33a02c' \
  --yAxisLabel 'Accessibility Signal' \
  --startLabel Start \
  --endLabel End \
  --plotHeight 12 \
  --plotWidth 13 \
  --plotFileFormat pdf \
  --dpi 720 \
  --plotTitle 'GV M2 CC Accessibility' \
  --perGroup &


# H3K4me3 ----------------------------------------------------------------------
WT_K4me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/WT/H3K4me3/deeptools/bamCoverage"
SCNT_K4me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/SCNT/H3K4me3/deeptools/bamCoverage" 

nohup computeMatrix scale-regions -p 20 -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
  -R "${DAR_dir}/DDARs.bed" "${DAR_dir}/Random_regions.bed" \
  -S "${WT_K4me3_bw}/MII_K4me3.bw" "${SCNT_K4me3_bw}/cc_K4me3.bw" \
  --samplesLabel WT_M2 SCNT_CC \
  --skipZeros \
  --outFileName M2_CC_K4me3_signal_around_DDARs.gz \
  --outFileNameMatrix M2_CC_K4me3_signal_around_DDARs.txt \
  --outFileSortedRegions M2_CC_K4me3_signal_around_DDARs.bed &

nohup plotProfile -m M2_CC_K4me3_signal_around_DDARs.gz \
  -out M2_CC_K4me3_signal_around_DDARs.pdf \
  --colors "#1f78b4" "#e31a1c" \
  --yAxisLabel 'H3K4me3 Signal' \
  --startLabel Start \
  --endLabel End \
  --plotHeight 12 \
  --plotWidth 13 \
  --plotFileFormat pdf \
  --dpi 720 \
  --plotTitle 'M2 vs CC H3K4me3' \
  --perGroup &


# H3K27me3 ---------------------------------------------------------------------
WT_K27me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/WT/H3K27me3/deeptools/bamCoverage"
SCNT_K27me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/SCNT/H3K27me3/deeptools/bamCoverage"   
ExE75_K27me3_bw="/home4/ssyi/Mouse_Placenta/Public_data/WMZ_E75EXE_H3K27me3NChIP/6.deeptools/bamCoverage/merge"

# CC/M2
nohup computeMatrix scale-regions -p 20 -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
  -R "${DAR_dir}/DDARs.bed" "${DAR_dir}/Random_regions.bed" \
  -S "${WT_K27me3_bw}/MII_K27me3.bw" "${SCNT_K27me3_bw}/cc_K27me3.bw" \
  --samplesLabel WT_M2 SCNT_CC \
  --skipZeros \
  --outFileName M2_CC_K27me3_signal_around_DDARs.gz \
  --outFileNameMatrix M2_CC_K27me3_signal_around_DDARs.txt \
  --outFileSortedRegions M2_CC_K27me3_signal_around_DDARs.bed &

nohup plotProfile -m M2_CC_K27me3_signal_around_DDARs.gz \
  -out M2_CC_K27me3_signal_around_DDARs.pdf \
  --colors "#1f78b4" "#e31a1c" \
  --yAxisLabel 'H3K27me3 Signal' \
  --startLabel Start \
  --endLabel End \
  --plotHeight 12 \
  --plotWidth 13 \
  --plotFileFormat pdf \
  --dpi 720 \
  --plotTitle 'M2 vs CC H3K27me3' \
  --perGroup &

# ExE7.5
nohup computeMatrix scale-regions -p 20 -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
  -R "${DAR_dir}/DDARs.bed" "${DAR_dir}/Random_regions.bed" \
  -S "${ExE75_K27me3_bw}/IVF_E75ExE.bw" "${ExE75_K27me3_bw}/SCNT_E75ExE.bw" \
  --samplesLabel IVF_E75ExE SCNT_E75ExE \
  --skipZeros \
  --outFileName E75ExE_K27me3_signal_around_DDARs.gz \
  --outFileNameMatrix E75ExE_K27me3_signal_around_DDARs.txt \
  --outFileSortedRegions E75ExE_K27me3_signal_around_DDARs.bed &

nohup plotProfile -m E75ExE_K27me3_signal_around_DDARs.gz \
  -out E75ExE_K27me3_signal_around_DDARs.pdf \
  --colors "#1f78b4" "#e31a1c" \
  --yAxisLabel 'H3K27me3 Signal' \
  --startLabel Start \
  --endLabel End \
  --plotHeight 12 \
  --plotWidth 13 \
  --plotFileFormat pdf \
  --dpi 720 \
  --plotTitle 'E75ExE WT vs SCNT H3K27me3' \
  --perGroup &


# H3K9me3 ----------------------------------------------------------------------
WT_K9me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/WT/H3K9me3/deeptools/bamCoverage"
SCNT_K9me3_bw="/home4/ssyi/Mouse_Placenta/codes/12.Modifications/Data/SCNT/H3K9me3/deeptools/bamCoverage"

nohup computeMatrix scale-regions -p 20 -b 3000 -a 3000 --regionBodyLength 6000 -bs 50 \
  -R "${DAR_dir}/DDARs.bed" "${DAR_dir}/Random_regions.bed" \
  -S "${WT_K9me3_bw}/Oocyte_K9me3.bw" "${SCNT_K9me3_bw}/CC_K9me3.bw" \
  --samplesLabel WT_M2 SCNT_CC \
  --skipZeros \
  --outFileName M2_CC_K9me3_signal_around_DDARs.gz \
  --outFileNameMatrix M2_CC_K9me3_signal_around_DDARs.txt \
  --outFileSortedRegions M2_CC_K9me3_signal_around_DDARs.bed &

nohup plotProfile -m M2_CC_K9me3_signal_around_DDARs.gz \
  -out M2_CC_K9me3_signal_around_DDARs.pdf \
  --colors "#1f78b4" "#e31a1c" \
  --yAxisLabel 'H3K9me3 Signal' \
  --startLabel Start \
  --endLabel End \
  --plotHeight 12 \
  --plotWidth 13 \
  --plotFileFormat pdf \
  --dpi 720 \
  --plotTitle 'M2 vs CC H3K9me3' \
  --perGroup &
