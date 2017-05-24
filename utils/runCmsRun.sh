#!/bin/bash -e

#$ -V
#$ -cwd
#$ -l os=sld6
#$ -l site=hh
#$ -P unihh2
#$ -m eas
#$ -M robin.aggleton@desy.de
#$ -t 1-17
#$ -l h_vmem=3G
#$ -l h_fsize=2G
#$ -l h_rt=8:00:00

# Run on BIRD: qsub runCmsRun.sh
# Don't forget to adjust the number of jobs after -t, and the walltime (h_rt)
# Don't lower vmem - will segfault otherwise

echo "Run job ${SGE_TASK_ID}"
cmsRun wrapper.py


