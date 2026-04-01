#!/bin/bash
#SBATCH -J gbm_download
#SBATCH -A FRANKELL-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --array=1-75
#SBATCH -o logs/%A_%a_dump_bam.log
#SBATCH --mail-type=ALL

module load samtools

# get srr id
srr=$( cat srr.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

echo "extracting bam for ${srr}"

apptainer exec ~/rds/rds-early-cancer_ev2-LH0AvU65IRI/_ENV_v2/singularity/sra-tools.sif sam-dump --ngc prj_41349.ngc ${srr}/${srr}.sra | samtools view -bS - > bams/${srr}.bam

echo "finished for ${srr}"
