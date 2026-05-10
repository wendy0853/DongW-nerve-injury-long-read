#!/bin/bash

# =============================================================================
# Run Cell Ranger count for single-nucleus RNA-seq alignment
# =============================================================================
#
# Purpose:
#   This script runs Cell Ranger count using the custom Mbp-split reference
#   generated with:
#
#       run_cellranger_mkref.sh
#
# Features enabled:
#   - BAM generation
#   - intronic read inclusion
#
# Docker image used:
#   jinlab/velocytoxcellranger:vs0.17.17x9.0.1
#
# Notes:
#   This script was originally executed on the WashU RIS HPC cluster using
#   Docker + LSF. Institution-specific submission commands were removed and
#   replaced with a generalized Cell Ranger command.
#
# =============================================================================

set -euo pipefail

# -----------------------------
# User-defined sample information
# -----------------------------

SAMPLE_ID="C0_1"                                            # <-- MODIFY HERE
FASTQ_SAMPLE_NAME="WT_Cell_Sciatic_1"                            # <-- MODIFY HERE

# -----------------------------
# User-defined paths
# -----------------------------

TRANSCRIPTOME_REF="/path/to/GRCm39_MbpSplit"                    # <-- MODIFY HERE
FASTQ_DIR="/path/to/FASTQ_directory"                            # <-- MODIFY HERE
OUTPUT_DIR="/path/to/output_directory"                          # <-- MODIFY HERE

# -----------------------------
# Resource settings
# -----------------------------

LOCAL_CORES=16                                                  # <-- MODIFY HERE IF NEEDED
LOCAL_MEM=120                                                   # <-- MODIFY HERE IF NEEDED

# -----------------------------
# Run Cell Ranger count
# -----------------------------

cellranger count \
  --id="${SAMPLE_ID}" \
  --transcriptome="${TRANSCRIPTOME_REF}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${FASTQ_SAMPLE_NAME}" \
  --output-dir="${OUTPUT_DIR}" \
  --create-bam=true \
  --include-introns=true \
  --localcores="${LOCAL_CORES}" \
  --localmem="${LOCAL_MEM}"
