##RUN CODE: bash bsub_sqanti3.sh {SAMPLE_NAME} {DIR} {SRB}

export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810/Active/References:/storage1/fs1/jin810/Active/References /storage1/fs1/jin810/Active/Projects:/storage1/fs1/jin810/Active/Projects /home/d.wendy:/home/d.wendy'
#-o %J.SQANTI3_C7_Injured_Sciatic_1.out -e %J.SQANTI_C7_Injured_Sciatic_1.err 
#C7_Injured_Sciatic_1
SAMPLE_NAME=$1
#C7_Injury
DIR=$2
#fastq_C7_Injured_1.fofn
SRB=$3
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub  -G compute-jin810-t3 -q subscription -sla jin810_t3 -n 16  -R 'rusage[mem=100GB]' -a 'docker(jinlab/sqanti3:vs1)' /bin/bash run_sqanti3.sh ${SAMPLE_NAME} ${DIR} ${SRB}

