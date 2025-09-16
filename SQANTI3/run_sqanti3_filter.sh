##!/bin/bash
#C7_Injured_Sciatic_1
SAMPLE_NAME=$1
#C7_Injury
DIR=$2

CLASS_TXT="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_classification.txt)"
CORRECTED_FASTA="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_corrected.fasta)"
CORRECTED_GTF="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_corrected.gtf)"
OUTPUT_DIR="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/filtered)"

# Create output directory
mkdir -p ${OUTPUT_DIR}
# Run SQANTI3 Filter
echo ${SAMPLE_NAME} ${CLASS_TXT} ${OUTPUT_DIR}

sqanti3_filter.py --sqanti_class ${CLASS_TXT} --filter_gtf ${CORRECTED_GTF} --dir ${OUTPUT_DIR} -o ${SAMPLE_NAME}_filtered
