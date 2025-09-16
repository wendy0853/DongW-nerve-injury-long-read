#!/bin/bash

#C7_Injured_Sciatic_1
SAMPLE_NAME=$1
#C7_Injury
DIR=$2
#fastq_C7_Injured_1.fofn
SRB=$3

ISOQUANT_GTF="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/IsoQuant_Results_updated/${SAMPLE_NAME}/${SAMPLE_NAME}.transcript_models.gtf)"
SHORT_READ_BAM="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/${SRB})"
OUTPUT_DIR="$(echo /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/${DIR}/SQANTI3_Output_updated/${SAMPLE_NAME})"


REFERENCE_GTF="/storage1/fs1/jin810/Active/References/GRCm39/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf"
REFERENCE_FASTA="/storage1/fs1/jin810/Active/References/GRCm39/mm38.fa"
POLYA_MOTIF_LIST="/opt2/sqanti3/5.3.6/SQANTI3-5.3.6/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
CAGE_PEAK_FILE="/opt2/sqanti3/5.3.6/SQANTI3-5.3.6/data/ref_TSS_annotation/mouse.refTSS_v3.1.GRCm39.bed"

# Create output directory
mkdir -p ${OUTPUT_DIR}

#pip install --no-cache-dir psutil

# Run SQANTI3 QC with CAGE peaks
#sqanti3_qc.py ${ISOQUANT_GTF} ${REFERENCE_GTF} ${REFERENCE_FASTA}   --SR_bam ${SHORT_READ_BAM}   --polyA_motif_list ${POLYA_MOTIF_LIST}   --CAGE_peak ${CAGE_PEAK_FILE}   -d ${OUTPUT_DIR}   --cpus 8 -o ${NAME}
echo ${SAMPLE_NAME} ${ISOQUANT_GTF} ${OUTPUT_DIR}
sqanti3_qc.py ${ISOQUANT_GTF} ${REFERENCE_GTF} ${REFERENCE_FASTA} --short_reads ${SHORT_READ_BAM}   --polyA_motif_list ${POLYA_MOTIF_LIST}   --CAGE_peak ${CAGE_PEAK_FILE}   -d ${OUTPUT_DIR}   --cpus 8 -o ${SAMPLE_NAME} 
