export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810/Active/References:/storage1/fs1/jin810/Active/References /storage1/fs1/jin810/Active/Projects:/storage1/fs1/jin810/Active/Projects /home/d.wendy:/home/d.wendy'

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub  -G compute-jin810-t3 -q subscription -sla jin810_t3 -n 16  -R 'rusage[mem=100GB]' -a 'docker(jinlab/sqanti3:vs1)' /bin/bash run_STAR_Alignment.sh
