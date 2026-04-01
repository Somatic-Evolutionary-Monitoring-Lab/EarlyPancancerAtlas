#!/bin/bash
#SBATCH -J gbm_download
#SBATCH -A FRANKELL-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH --array=1-75
#SBATCH -o logs/%A_%a_get_bam.log
#SBATCH --mail-type=ALL

# get srr id
srr=$( cat srr.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

apptainer exec ~/rds/rds-early-cancer_ev2-LH0AvU65IRI/_ENV_v2/singularity/sra-tools.sif prefetch -p --ngc prj_41349.ngc $srr
