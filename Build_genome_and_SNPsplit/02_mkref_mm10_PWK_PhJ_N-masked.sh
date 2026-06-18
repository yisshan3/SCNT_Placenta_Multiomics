#!/bin/sh

# ==============================================================================
# Description: Custom reference building for Cell Ranger ARC using the PWK/PhJ N-masked genome. Adapts FASTA headers and filters GTF biotypes.
# ==============================================================================

cd /home1/ssyi/annotate/mm/SNP/cellranger-arc

# Version
# cellranger-arc-2.0.2

# Genome metadata
genome="mm10_PWK_PhJ_N-masked"
version="release98"


# Set up source and build directories
build="${genome}-${version}-build"
mkdir -p "$build"


# Download source files if they do not exist in the reference_sources folder
source="${genome}-${version}-reference-sources"
mkdir -p "$source"


fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
fasta_in="${source}/mm10_PWK_PhJ_N-masked.fasta"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.vM23.primary_assembly.annotation.gtf"
motifs_url="https://jaspar.genereg.net/download/data/2018/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
motifs_in="${source}/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"


# Check for missing input files and download if necessary
if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi
if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi


# Modify sequence headers in the Ensembl FASTA to match the GENCODE format.
# Unplaced and unlocalized sequences have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"

# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" prefix to autosomes and sex chromosomes
# 3. Handle the mitochondrial chromosome (MT -> chrM)
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages.
#
# Input GTF:
#     ... gene_id "ENSMUSG00000102693.1"; ...
# Output GTF:
#     ... gene_id "ENSMUSG00000102693"; gene_version "1"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"

# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
chr=\
"(chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|\
chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|\
chrX|chrY|chrM)"

cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    | grep -E "$chr" \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lncRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"

GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""


# Construct the gene ID allowlist. 
# We filter the list of all transcripts based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"

# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"

# Filter to the gene allowlist
# Note: The modified GTF now only includes standard chromosomes (chr1-19, X, Y, M) 
# and the filtered gene lists.
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"


# Change motif headers so the human-readable motif name precedes the motif
# identifier. E.g., ">MA0004.1 Arnt" -> ">Arnt_MA0004.1".
motifs_modified="$build/$(basename "$motifs_in").modified"
awk '{
    if ( substr($1, 1, 1) == ">" ) {
        print ">" $2 "_" substr($1,2)
    } else {
        print
    }
}' "$motifs_in" > "$motifs_modified"


# Create a config file for mkref
config_in="${build}/config"
echo """{
    organism: \"Mus_musculus\"
    genome: [\""$genome"\"]
    input_fasta: [\""$fasta_modified"\"]
    input_gtf: [\""$gtf_filtered"\"]
    input_motifs: \""$motifs_modified"\"
    non_nuclear_contigs: [\"chrM\"]
}""" > "$config_in"


# Create reference package
cellranger-arc mkref --ref-version="$version" \
    --config="$config_in"
    