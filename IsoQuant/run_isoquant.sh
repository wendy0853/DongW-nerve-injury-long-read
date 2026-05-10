#!/bin/bash

# =========================
# IsoQuant Long-Read Pipeline
# =========================

# Activate conda environment
export PATH=/miniconda/bin/:$PATH
source /miniconda/etc/profile.d/conda.sh
conda activate isoquant

# =========================
# Input Parameters
# =========================

SAMPLE_NAME="C0_Sciatic_1" #<-- Change for your sample

FASTQ_DIR="/path/to/fastq_files"
OUTPUT_DIR="/path/to/output_directory"

REFERENCE_FASTA="/path/to/GRCm39.genome.fa.gz"
REFERENCE_GTF="/path/to/gencode.vM38.annotation.gtf.gz"

# =========================
# Merge FASTQ Files
# =========================

cat ${FASTQ_DIR}/*.fastq | gzip -c > ${OUTPUT_DIR}/${SAMPLE_NAME}.fastq.gz

# =========================
# Run IsoQuant
# =========================

isoquant.py \
-d nanopore \
--stranded forward \
--fastq ${OUTPUT_DIR}/${SAMPLE_NAME}.fastq.gz \
--reference ${REFERENCE_FASTA} \
--genedb ${REFERENCE_GTF} \
--complete_genedb \
--sqanti_output \
--check_canonical \
--count_exons \
--report_novel_unspliced true \
--output ${OUTPUT_DIR} \
--prefix ${SAMPLE_NAME}
