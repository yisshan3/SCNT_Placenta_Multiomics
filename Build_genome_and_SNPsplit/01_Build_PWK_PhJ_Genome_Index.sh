#!/bin/sh

# ==============================================================================
# Description: Download reference genome/annotation, generate PWK/PhJ N-masked genome, and build indices for Cell Ranger ARC, HISAT2, and Bowtie2.
# ==============================================================================

work_dir="/home1/ssyi/annotate/mm/SNP"
cd "${work_dir}" || exit


# Download Reference Genome (mm10/GRCm38) --------------------------------------
mkdir -p mm10_fasta
cd mm10_fasta || exit

# Download FASTA file
wget -c http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/README
gunzip -v Mus_musculus.GRCm38.dna.primary_assembly.fa.gz


# Download Genome Annotation (Gencode vM23) ------------------------------------
cd "${work_dir}"
mkdir -p mm10_gtf
cd mm10_gtf || exit

# Download GTF file
wget -c http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
wget -c http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/MD5SUMS
wget -c http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/_README.TXT
# Verify MD5 checksum
md5sum -c MD5SUMS
gunzip -v gencode.vM23.primary_assembly.annotation.gtf.gz


# PWK/PhJ N-masked Genome Preparation ------------------------------------------
cd "${work_dir}"
# Download PWK/PhJ VCF files
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz.md5 &

mkdir -p PWK_PhJ_Single_strain
cd PWK_PhJ_Single_strain || exit

# Run SNPsplit genome preparation to generate the N-masked genome
nohup SNPsplit_genome_preparation --vcf_file ../PWK_PhJ.mgp.v5.snps.dbSNP142.vcf.gz \
    --reference_genome "${work_dir}/mm10_fasta/" \
    --strain PWK_PhJ --nmasking --genome_build GRCm38 &
# Note: Generates all_SNPs_PWK_PhJ_GRCm38.txt.gz, reports, and N-masked fasta files

gzip -d all_SNPs_PWK_PhJ_GRCm38.txt.gz

# Add "chr" prefix to chromosome names in the SNP list
awk 'BEGIN{FS=OFS="\t"}$2="chr"$2{print $0}' all_SNPs_PWK_PhJ_GRCm38.txt > all_SNPs_PWK_PhJ_GRCm38_modified.txt

# Concatenate standard chromosomes (chr1-19, X, Y, M)
cd PWK_PhJ_N-masked
cat *.N-masked.fa > mm10_PWK_PhJ_N-masked.fasta   


# Build Cell Ranger ARC Reference ----------------------------------------------
cd "${work_dir}"
mkdir -p cellranger-arc/mm10_PWK_N-masked
cd cellranger-arc/mm10_PWK_N-masked || exit

# Set up reference sources
mkdir -p mm10_PWK_PhJ_N-masked-release98-build mm10_PWK_PhJ_N-masked-release98-reference-sources
cd mm10_PWK_PhJ_N-masked-release98-reference-sources || exit

# Download JASPAR motifs
nohup wget -c --no-check-certificate https://jaspar.genereg.net/download/data/2018/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt &
ln -s "${work_dir}/mm10_gtf/gencode.vM23.primary_assembly.annotation.gtf" ./
ln -s "${work_dir}/PWK_PhJ_Single_strain/PWK_PhJ_N-masked/mm10_PWK_PhJ_N-masked.fasta" ./

# Execute 02_mkref_mm10_PWK_PhJ_N-masked.sh for Cell Ranger ARC reference
nohup sh 02_mkref_mm10_PWK_PhJ_N-masked.sh &


# Build HISAT2 Index for N-masked Genome ---------------------------------------
cd "${work_dir}"
mkdir -p PWK_PhJ_N-masked_hisat2_index
cd PWK_PhJ_N-masked_hisat2_index || exit

nohup hisat2_extract_exons.py "${work_dir}/PWK_PhJ_Single_strain/gencode.vM23.primary_assembly.annotation.gtf.filtered" > exons.gtf &
nohup hisat2_extract_splice_sites.py "${work_dir}/PWK_PhJ_Single_strain/gencode.vM23.primary_assembly.annotation.gtf.filtered" > splice_sites.gtf &
nohup hisat2-build -p 10 --ss splice_sites.gtf --exon exons.gtf "${work_dir}/PWK_PhJ_Single_strain/mm10_PWK_PhJ_N-masked.fasta.modified" PWK_PhJ_N-masked &


# Build Bowtie2 Index for N-masked Genome --------------------------------------
cd "${work_dir}"
mkdir -p PWK_PhJ_N-masked_bowtie2_index
cd PWK_PhJ_N-masked_bowtie2_index || exit

nohup bowtie2-build -f "${work_dir}/PWK_PhJ_Single_strain/mm10_PWK_PhJ_N-masked.fasta.modified" --threads 15 PWK_PhJ_N-masked &

