export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810:/storage1/fs1/jin810 /home/d.wendy:/home/d.wendy'

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub \
  -G compute-jin810-t3 \
  -q subscription \
  -sla jin810_t3 \
  -n 16 \
  -M 200GB \
  -R 'span[hosts=1] rusage[mem=150GB]' \
  -a 'docker(jinlab/velocytoxcellranger:vs0.17.17x9.0.1)' \
  cellranger mkref \
    --genome=GRCm39_MbpSplit \
    --fasta=/storage1/fs1/jin810/Active/References/GRCm39_10X/refdata-gex-GRCm39-2024-A/fasta/genome.fa \
    --genes=/storage1/fs1/jin810/Active/References/GRCm39_10X/genes.MbpSplit.gtf
