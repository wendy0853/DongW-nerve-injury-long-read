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

sqanti3_filter.py ml --sqanti_class ${CLASS_TXT} --filter_gtf ${CORRECTED_GTF} --dir ${OUTPUT_DIR} -o ${SAMPLE_NAME}_filtered

sqanti3_rescue.py rules \
    --filter_class example/rules_filter_output/UHR_chr22_RulesFilter_result_classification.txt \
    --refGTF data/reference/gencode.v38.basic_chr22.gtf \
    -refFasta data/reference/GRCh38.p13_chr22.fasta \
    --rescue_isoforms example/SQANTI3_QC_output/UHR_chr22_corrected.fasta \
    --rescue_gtf example/rules_filter_output/UHR_chr22.filtered.gtf \
    --refClassif data/reference/gencode.v38.basic_chr22_classification.txt \
    --mode full \
    --json_filter src/utilities/filter/filter_default.json \
    --dir example/rescue_full_rules --output UHR_chr22 
