# Start interactive session
bsub -Is -G compute-jin810-t3 -q subscription -sla jin810_t3 -n 8 \
    -R 'gpuhost rusage[mem=64GB]' \
    -gpu 'num=1:j_exclusive=yes' \
    -a 'docker(biocontainers/samtools:v1.9-4-deb_cv1)' \
/bin/bash


# 1) extract Jun region (+/- 5kb padding)
samtools view -b \
  C7_Injured_Sciatic_3_C7_Injured_Sciatic_3_98f996_7ac08f_ae73b0.bam \
  chr4:94932271-94945459 \
  > C7_Injured_Sciatic_3_Jun.bam

# 2) index the subset BAM
samtools index C7_Injured_Sciatic_3_Jun.bam


