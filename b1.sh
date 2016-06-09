#!/bin/bash
#
#SBATCH --job-name=shelx

cd /mnt/das-gpfs/work/p11206/tubu_2.5A/shelx_v2c100k 

hostname
echo Start running
date
python ../shelx_batch.py ../24turns.ahkl $1 $2 $3 
date

