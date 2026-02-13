REF_DIR="/storage1/fs1/jin810/Active/References/GRCm39"
FA="${REF_DIR}/GRCm39.genome.fa"
GTF="${REF_DIR}/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf"
GENOME_DIR="${REF_DIR}/STAR_gencode_vM38_GRCm39"

THREADS="${THREADS:-16}"
SJDB_OVERHANG="${SJDB_OVERHANG:-149}" 

mkdir -p "${GENOME_DIR}"

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "${REF_DIR}:/ref:ro" \
  -v "${GENOME_DIR}:/genome" \
  "${STAR_IMG}" \
  STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir /genome \
    --genomeFastaFiles "/ref/$(basename "${FA}")" \
    --sjdbGTFfile "/ref/$(basename "${GTF}")" \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbOverhang "${SJDB_OVERHANG}"

STAR --runThreadN 16 \
  --genomeDir star_newref \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix sample.newref.
