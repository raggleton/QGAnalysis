#!/bin/bash -e

#$ -V
#$ -cwd
#$ -l os=sld6
#$ -l site=hh
#$ -P unihh2
#$ -m eas
#$ -M aggleton@desy.de
#$ -t 1-10
#$ -l h_vmem=2G
#$ -l h_fsize=2G
#$ -l h_rt=3:00:00

echo "Run job ${SGE_TASK_ID}"
cmsRun wrapper.py


