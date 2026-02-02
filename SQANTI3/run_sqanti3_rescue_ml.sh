##!/bin/bash
#C7_Injured_Sciatic_1
SAMPLE_NAME=$1
#C7_Injury
DIR=$2

CLASS_TXT="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_classification.txt)"
CORRECTED_FASTA="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_corrected.fasta)"
CORRECTED_GTF="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/${SAMPLE_NAME}_corrected.gtf)"
OUTPUT_DIR="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME}/filtered_ml)"
REFERENCE_GTF="/storage1/fs1/jin810/Active/References/GRCm39/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"
REFERENCE_FASTA="/storage1/fs1/jin810/Active/References/GRCm39/GRCm39.genome.fa"
POLYA_MOTIF_LIST="/opt2/sqanti3/5.3.6/SQANTI3-5.3.6/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
CAGE_PEAK_FILE="/opt2/sqanti3/5.3.6/SQANTI3-5.3.6/data/ref_TSS_annotation/mouse.refTSS_v3.1.GRCm39.bed"

# Create output directory
mkdir -p ${OUTPUT_DIR}
# Run SQANTI3 Filter
echo ${SAMPLE_NAME} ${CLASS_TXT} ${OUTPUT_DIR}

sqanti3_rescue.py ml --filter_class  ${OUTPUT_DIR} \
                        --refGTF ${REFERENCE_GTF} \
                        --refFasta ${REFERENCE_FASTA} \
                        --mode automatic \
                        --dir ${OUTPUT_DIR}/rescue_automatic --output ${SAMPLE_NAME}_rescue
