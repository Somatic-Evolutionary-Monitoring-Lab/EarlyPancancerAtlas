#!/bin/bash
#SBATCH -J gbm_download
#SBATCH -A FRANKELL-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --array=1-150
#SBATCH -o logs/%A_%a_dump_bam.log
#SBATCH --mail-type=ALL

module load fastqc

# get srr id
fq=$( cat fastqs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#  fastqc
fastqc --threads 4 $fq

