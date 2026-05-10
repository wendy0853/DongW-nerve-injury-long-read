#!/bin/bash

# =============================================================================
# Run CellBender ambient RNA correction
# =============================================================================
#
# Purpose:
#   This script runs CellBender remove-background to remove ambient RNA
#   contamination from single-nucleus RNA-seq count matrices.
#
# Docker image used:
#   us.gcr.io/broad-dsde-methods/cellbender:0.3.0
#
# Notes:
#   This workflow was originally run on the WashU RIS HPC cluster using
#   Docker + GPU-enabled LSF submission. Institution-specific submission
#   commands were removed and replaced with a generalized CellBender command.
#
# Requirements:
#   - GPU-enabled environment recommended
#   - CUDA-compatible installation if not using Docker
#
# =============================================================================

set -euo pipefail

# -----------------------------
# User-defined input/output files
# -----------------------------

INPUT_H5="/path/to/raw_feature_bc_matrix.h5"               # <-- MODIFY HERE
OUTPUT_H5="/path/to/cellbender_output_file.h5"             # <-- MODIFY HERE

# -----------------------------
# User-adjustable parameters
# -----------------------------

EXPECTED_CELLS=30000                                       # <-- MODIFY HERE DEPENDING ON CELLS
LEARNING_RATE=0.000001                                     # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Run CellBender
# -----------------------------

cellbender remove-background \
  --cuda \
  --input "${INPUT_H5}" \
  --output "${OUTPUT_H5}" \
  --expected-cells "${EXPECTED_CELLS}" \
  --learning-rate "${LEARNING_RATE}"
