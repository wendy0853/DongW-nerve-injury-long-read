#!/bin/bash

# =========================
# SQANTI3 Rules-Based Filtering
# Executed inside Docker container:
# jinlab/sqanti3:vs1
# =========================

set -euo pipefail

# Usage:
# bash run_sqanti3_filter.sh SAMPLE_NAME

SAMPLE_NAME=$1

# =========================
# User-defined paths
# =========================

SQANTI3_DIR="/path/to/SQANTI3/output/${SAMPLE_NAME}"
OUTPUT_DIR="${SQANTI3_DIR}/filtered_rules"

CLASS_TXT="${SQANTI3_DIR}/${SAMPLE_NAME}_classification.txt"
CORRECTED_GTF="${SQANTI3_DIR}/${SAMPLE_NAME}_corrected.gtf"

# =========================
# Run SQANTI3 filter
# =========================

mkdir -p "${OUTPUT_DIR}"

sqanti3_filter.py \
rules \
--sqanti_class "${CLASS_TXT}" \
--filter_gtf "${CORRECTED_GTF}" \
--dir "${OUTPUT_DIR}" \
-o "${SAMPLE_NAME}_filtered"
