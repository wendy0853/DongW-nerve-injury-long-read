#!/bin/bash

# =========================
# SQANTI3 QC
# Executed inside Docker container:
# jinlab/sqanti3:vs1
# =========================

set -euo pipefail

# Usage:
# bash run_sqanti3_qc.sh SAMPLE_NAME

SAMPLE_NAME=$1

# =========================
# User-defined paths
# =========================

ISOQUANT_GTF="/path/to/IsoQuant/output/${SAMPLE_NAME}.transcript_models.gtf"
SHORT_READS="/path/to/short_read_alignments_or_fofn/${SAMPLE_NAME}.fofn"
OUTPUT_DIR="/path/to/SQANTI3/output/${SAMPLE_NAME}"

REFERENCE_GTF="/path/to/reference/gencode.vM38.annotation.gtf"
REFERENCE_FASTA="/path/to/reference/GRCm39.genome.fa"

POLYA_MOTIF_LIST="/path/to/SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
CAGE_PEAK_FILE="/path/to/SQANTI3/data/ref_TSS_annotation/mouse.refTSS_v3.1.GRCm39.bed"

# =========================
# Run SQANTI3 QC
# =========================

mkdir -p "${OUTPUT_DIR}"

sqanti3_qc.py \
"${ISOQUANT_GTF}" \
"${REFERENCE_GTF}" \
"${REFERENCE_FASTA}" \
--short_reads "${SHORT_READS}" \
--polyA_motif_list "${POLYA_MOTIF_LIST}" \
--CAGE_peak "${CAGE_PEAK_FILE}" \
-d "${OUTPUT_DIR}" \
--cpus 8 \
-o "${SAMPLE_NAME}"
