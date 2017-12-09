#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -o log
#$ -e log
#$ -t 1-9504
#$ -tc 700
#$ -l h_vmem=8G
#$ -l long

hostname
date
file="lm_coef_p/$SGE_TASK_ID.txt"
Rscript run_lm_peer_L1_noPerm.R data4lm_peer/$SGE_TASK_ID.txt
date

