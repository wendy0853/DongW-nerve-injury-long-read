#!/bin/bash

# =========================
# SQANTI3 Automatic Rescue
# Executed inside Docker container:
# jinlab/sqanti3:vs1
# =========================

set -euo pipefail

# Usage:
# bash run_sqanti3_rescue_rules_automatic.sh SAMPLE_NAME

SAMPLE_NAME=$1

# =========================
# User-defined paths
# =========================

SQANTI3_DIR="/path/to/SQANTI3/output/${SAMPLE_NAME}"
FILTER_DIR="${SQANTI3_DIR}/filtered_rules"
OUTPUT_DIR="${FILTER_DIR}/rescue_automatic"

REFERENCE_GTF="/path/to/reference/gencode.vM38.annotation.gtf"
REFERENCE_FASTA="/path/to/reference/GRCm39.genome.fa"

FILTERED_CLASS="${FILTER_DIR}/${SAMPLE_NAME}_filtered_RulesFilter_result_classification.txt"

# =========================
# Run SQANTI3 automatic rescue
# =========================

mkdir -p "${OUTPUT_DIR}"

sqanti3_rescue.py \
rules \
--filter_class "${FILTERED_CLASS}" \
--refGTF "${REFERENCE_GTF}" \
--refFasta "${REFERENCE_FASTA}" \
--mode automatic \
--dir "${OUTPUT_DIR}" \
--output "${SAMPLE_NAME}_rescue"
