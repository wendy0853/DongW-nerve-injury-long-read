export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /home/d.wendy:/home/d.wendy'

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id WT_Mbp_1 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20241009_GEMX_Julie/FASTQ \
--sample WT_Nuc_Sciatic_1 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/WT_Mbp_1 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id WT_Mbp_2 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20241009_GEMX_Julie/FASTQ \
--sample WT_Nuc_Sciatic_2 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/WT_Mbp_2 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id TH_Mbp_1 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20241009_GEMX_Julie/FASTQ \
--sample TH_Nuc_Sciatic_1 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/TH_Mbp_1 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id TH_Mbp_2 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20241009_GEMX_Julie/FASTQ \
--sample TH_Nuc_Sciatic_2 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/TH_Mbp_2 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id SOD1_Mbp_1 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20250925_SOD1_G93A_SN_NS0_NS1/FASTQ \
--sample GEMX_M_NS1_SOD1_G93A_Sciatic \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/SOD1_Mbp_1 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id SOD1_Mbp_2 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/Projects/snRNAseq/20251217_GEMX_SOD1_G93A_NS0_Fem_Sural_NS01_SN/FASTQ \
--sample LIB128344-DIL01-DIL01_23GHN2LT4 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260212_Golli_Compact_Mbp_Alignment/SOD1_Mbp_2 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-jin810-t3 -q subscription -sla jin810_t3 -M 200GB -R 'rusage[mem=150GB]' \
-a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
cellranger count \
--id C3_1 \
--transcriptome /storage1/fs1/jin810/Active/References/GRCm39_10X/GRCm39_MbpSplit \
--fastqs /storage1/fs1/jin810/Active/U19_Data_Core/1_Project/Milbrandt_Lab/Datasets/isnat_data_22Feb2026/SRR18355056 \
--sample Sample_1181-RG-01 \
--output-dir /storage1/fs1/jin810/Active/Projects/snRNAseq/20260223_Golli_Compact_Mbp_iSNAT/C3_1 \
--create-bam=true --include-introns=true --localcores=16 --localmem=120


LSF_DOCKER_ENTRYPOINT=/bin/sh bsub -Is -G compute-jin810 -q general-interactive -n 4 -a 'docker(elle/basic:vs5)' /bin/sh
