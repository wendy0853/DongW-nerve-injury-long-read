# 0) sanity: make sure the BAM is indexed
samtools index sample.sorted.bam

# 1) extract Jun region (+/- 5kb padding)
samtools view -b \
  sample.sorted.bam \
  chr4:94932271-94945459 \
  > sample_Jun.bam

# 2) index the subset BAM
samtools index sample_Jun.bam


bsub -Is -G compute-jin810 -q general-interactive -n 4 \
    -R 'gpuhost rusage[mem=64GB]' \
    -gpu 'num=1:j_exclusive=yes' \
    -a 'docker(us.gcr.io/broad-dsde-methods/cellbender:0.3.0)' \
/bin/bash
