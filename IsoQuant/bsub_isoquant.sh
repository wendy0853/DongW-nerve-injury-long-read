export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810/Active/References:/storage1/fs1/jin810/Active/References /storage1/fs1/jin810/Active/Projects:/storage1/fs1/jin810/Active/Projects /home/d.wendy:/home/d.wendy'

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q general -n 16 -M 64GB -R "rusage[mem=64GB] span[hosts=1]" -G compute-jin810 -a 'docker(etycksen/isoquant:latest)' /bin/bash run_isoquant.sh
