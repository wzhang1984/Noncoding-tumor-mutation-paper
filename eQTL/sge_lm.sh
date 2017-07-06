#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -o log
#$ -e log
#$ -t 1-5346
#$ -tc 700
#$ -l h_vmem=8G
#$ -l long

hostname
date
file="lm_coef_p/$SGE_TASK_ID.txt"
if [ -f "$file" ]
then
    echo "$file found."
else
    echo "$file not found."
    Rscript run_lm_peer_L1_noPerm.R data4lm_peer/$SGE_TASK_ID.txt
fi
date
