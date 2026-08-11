#!/bin/bash

# =============================================================================
# Create modified GRCm39 Cell Ranger reference annotation with split Mbp models
# =============================================================================
#
# Purpose:
#   This script modifies the 10x Genomics GRCm39 Cell Ranger reference GTF by:
#
#     1. Removing the original canonical Mbp gene annotation
#     2. Creating a custom Mbp_Golli gene model
#     3. Creating a custom Mbp_Classic gene model
#     4. Writing a modified GTF for downstream Cell Ranger mkref
#
# Input:
#   - refdata-gex-GRCm39-2024-A.tar.gz
#
# Output:
#   - genes.MbpSplit.gtf
#
# Notes:
#   This script was originally run on the WashU RIS environment but has
#   been generalized for reuse. Update paths marked with "<-- MODIFY HERE".
#
# =============================================================================

set -euo pipefail

# -----------------------------
# User-defined inputs
# -----------------------------

REFERENCE_TAR="/path/to/refdata-gex-GRCm39-2024-A.tar.gz"     # <-- MODIFY HERE
WORKDIR="/path/to/GRCm39_2024A_MbpSplitWork"                 # <-- MODIFY HERE

MBP_GENE_ID="ENSMUSG00000041607"                             # Mbp Ensembl gene ID

mkdir -p "${WORKDIR}"

# -----------------------------
# Extract Cell Ranger reference
# -----------------------------

tar -xzf "${REFERENCE_TAR}" -C "${WORKDIR}"

REFDIR="$(find "${WORKDIR}" -maxdepth 1 -type d -name 'refdata-gex-GRCm39-2024-A*' | head -n 1)"

if [[ -z "${REFDIR}" ]]; then
  echo "ERROR: Could not find extracted Cell Ranger reference directory."
  exit 1
fi

echo "Reference directory: ${REFDIR}"

GTF_GZ="${REFDIR}/genes/genes.gtf.gz"
GTF_IN="${REFDIR}/genes/genes.gtf"
GENOME_FASTA="${REFDIR}/fasta/genome.fa"

if [[ ! -f "${GTF_GZ}" ]]; then
  echo "ERROR: Missing input GTF: ${GTF_GZ}"
  exit 1
fi

if [[ ! -f "${GENOME_FASTA}" ]]; then
  echo "ERROR: Missing genome FASTA: ${GENOME_FASTA}"
  exit 1
fi

gunzip -c "${GTF_GZ}" > "${GTF_IN}"

echo "Uncompressed GTF written to: ${GTF_IN}"

# -----------------------------
# Confirm Mbp annotation exists
# -----------------------------

echo "Checking for Mbp annotation..."
grep -m 1 'gene_name "Mbp"' "${GTF_IN}" || {
  echo "ERROR: Could not find Mbp annotation in ${GTF_IN}"
  exit 1
}

# -----------------------------
# Define exon IDs used for custom Mbp models
# -----------------------------
#
# These exon IDs were selected to distinguish Golli-associated and classic
# myelin-associated Mbp transcript models in the custom Cell Ranger reference.

GOLLI_PATTERNS="${WORKDIR}/MbpGolli.patterns.txt"
CLASSIC_PATTERNS="${WORKDIR}/MbpClassic.patterns.txt"

cat > "${GOLLI_PATTERNS}" << 'EOF'
exon_id "ENSMUSE00000571789"
exon_id "ENSMUSE00001256352"
exon_id "ENSMUSE00001241294"
EOF

cat > "${CLASSIC_PATTERNS}" << 'EOF'
exon_id "ENSMUSE00000376561"
exon_id "ENSMUSE00001457516"
exon_id "ENSMUSE00000349836"
EOF

# -----------------------------
# Verify selected exon IDs exist
# -----------------------------

echo "Golli exon hits:"
grep -F -f "${GOLLI_PATTERNS}" "${GTF_IN}" | head -n 5 || true
echo "Number of Golli exon hits:"
grep -F -f "${GOLLI_PATTERNS}" "${GTF_IN}" | wc -l

echo "Classic exon hits:"
grep -F -f "${CLASSIC_PATTERNS}" "${GTF_IN}" | head -n 5 || true
echo "Number of Classic exon hits:"
grep -F -f "${CLASSIC_PATTERNS}" "${GTF_IN}" | wc -l

# -----------------------------
# Build minimal exon GTFs
# -----------------------------

GOLLI_EXONS="${WORKDIR}/Mbp_Golli.exons.min.gtf"
CLASSIC_EXONS="${WORKDIR}/Mbp_Classic.exons.min.gtf"

grep -F -f "${GOLLI_PATTERNS}" "${GTF_IN}" |
  awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
    ex="NA";
    if (match($9, /exon_id "([^"]+)"/, a)) ex=a[1];
    key=$1 FS $4 FS $5 FS $7;
    if (!(key in seen)) {
      seen[key]=1;
      print $1,"custom","exon",$4,$5,".",$7,".",
        "gene_id \"Mbp_Golli\"; transcript_id \"Mbp_Golli_tx\"; gene_name \"Mbp_Golli\"; transcript_name \"Mbp_Golli_tx\"; exon_id \""ex"\";"
    }
  }' > "${GOLLI_EXONS}"

grep -F -f "${CLASSIC_PATTERNS}" "${GTF_IN}" |
  awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {
    ex="NA";
    if (match($9, /exon_id "([^"]+)"/, a)) ex=a[1];
    key=$1 FS $4 FS $5 FS $7;
    if (!(key in seen)) {
      seen[key]=1;
      print $1,"custom","exon",$4,$5,".",$7,".",
        "gene_id \"Mbp_Classic\"; transcript_id \"Mbp_Classic_tx\"; gene_name \"Mbp_Classic\"; transcript_name \"Mbp_Classic_tx\"; exon_id \""ex"\";"
    }
  }' > "${CLASSIC_EXONS}"

# -----------------------------
# Validate custom exon GTFs
# -----------------------------

awk -F'\t' '$0 !~ /^#/ && NF!=9{print "BAD",FNR,NF; exit 1}' "${GOLLI_EXONS}"
awk -F'\t' '$0 !~ /^#/ && NF!=9{print "BAD",FNR,NF; exit 1}' "${CLASSIC_EXONS}"

echo "Unique genomic exon intervals retained:"
echo -n "  Golli:   "
wc -l "${GOLLI_EXONS}"
echo -n "  Classic: "
wc -l "${CLASSIC_EXONS}"

# -----------------------------
# Build gene and transcript records
# -----------------------------

make_gene_tx_simple () {
  local exongtf="$1"
  local gid="$2"
  local gname="$3"
  local txid="$4"
  local txname="$5"

  local chrom
  local strand
  local start
  local end

  chrom=$(awk -F'\t' 'NR==1{print $1}' "${exongtf}")
  strand=$(awk -F'\t' 'NR==1{print $7}' "${exongtf}")
  start=$(awk -F'\t' 'BEGIN{m=1e18} {if($4<m)m=$4} END{print m}' "${exongtf}")
  end=$(awk -F'\t' 'BEGIN{m=0} {if($5>m)m=$5} END{print m}' "${exongtf}")

  printf "%s\tcustom\tgene\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; gene_type \"protein_coding\";\n" \
    "${chrom}" "${start}" "${end}" "${strand}" "${gid}" "${gname}"

  printf "%s\tcustom\ttranscript\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n" \
    "${chrom}" "${start}" "${end}" "${strand}" "${gid}" "${txid}" "${gname}" "${txname}"
}

GOLLI_GENE_TX="${WORKDIR}/Mbp_Golli.gene_tx.gtf"
CLASSIC_GENE_TX="${WORKDIR}/Mbp_Classic.gene_tx.gtf"

make_gene_tx_simple "${GOLLI_EXONS}" \
  "Mbp_Golli" "Mbp_Golli" "Mbp_Golli_tx" "Mbp_Golli_tx" \
  > "${GOLLI_GENE_TX}"

make_gene_tx_simple "${CLASSIC_EXONS}" \
  "Mbp_Classic" "Mbp_Classic" "Mbp_Classic_tx" "Mbp_Classic_tx" \
  > "${CLASSIC_GENE_TX}"

# -----------------------------
# Remove original Mbp and append custom split annotations
# -----------------------------

GENES_NO_MBP="${WORKDIR}/genes.noMbp.gtf"
GENES_MBP_SPLIT="${WORKDIR}/genes.MbpSplit.gtf"

awk -F'\t' -v gid="${MBP_GENE_ID}" '
  $0 ~ /^#/ {print; next}
  $9 !~ ("gene_id \""gid"\"") {print}
' "${GTF_IN}" > "${GENES_NO_MBP}"

cat "${GENES_NO_MBP}" \
    "${GOLLI_GENE_TX}" "${GOLLI_EXONS}" \
    "${CLASSIC_GENE_TX}" "${CLASSIC_EXONS}" \
  > "${GENES_MBP_SPLIT}"

# -----------------------------
# Structural validation
# -----------------------------

awk -F'\t' '
  $0 ~ /^#/ {next}
  NF!=9 {print "BAD", NR, NF; bad=1}
  END {if (bad) exit 1; else print "OK: all non-comment lines have 9 columns"}
' "${GENES_MBP_SPLIT}"

echo "Modified GTF written to: ${GENES_MBP_SPLIT}"
echo "Genome FASTA for mkref: ${GENOME_FASTA}"

