export PATH=/miniconda/bin/:$PATH && \
source /miniconda/etc/profile.d/conda.sh && \
conda activate isoquant && \
cat /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/C0_Control/Long_read_raw/fastq_Dorado/C0_Sciatic_1/*.fastq | gzip -c > ./storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/C0_Control/Long_read_raw/fastq_Dorado/C0_Sciatic_1/C0_Sciatic_1.fastq.gz && \
isoquant.py -d nanopore --stranded forward \
--fastq /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/C0_Control/Long_read_raw/fastq_Dorado/C0_Sciatic_1/C0_Sciatic_1.fastq.gz \
--reference /storage1/fs1/jin810/Active/References/GRCm39/GRCm39.genome.fa.gz \
--genedb /storage1/fs1/jin810/Active/References/GRCm39/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz \
--complete_genedb \
--sqanti_output --check_canonical --count_exons \
--report_novel_unspliced true \
--output /storage1/fs1/jin810/Active/Projects/RNAseq/Long_read_RNAseq_Crush_Injury/Long_read/C0_Control/IsoQuant_Results_updated \
--prefix C0_Sciatic_1
