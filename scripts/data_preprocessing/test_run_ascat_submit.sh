#!/bin/sh

#SBATCH -J run_alignment
#SBATCH -A FRANKELL-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --array=1-68
#SBATCH -o logs/ascat_test_%A_%a.log
#SBATCH --mail-type=ALL

cram=$( cat samples_to_run.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
sample=$( echo $cram | cut -d. -f1 )

conda activate ~/rds/rds-early-cancer_ev2-LH0AvU65IRI/_ENV_v2/conda/ascat_hts
export PATH=~/rds/rds-early-cancer_ev2-LH0AvU65IRI/_ENV_v2/conda/ascat_hts/bin:$PATH

# run script
Rscript ascat_germline_detect.R --bam $cram --sample_id $sample --outdir /home/kh723/rds/rds-early-cancer_ev2-LH0AvU65IRI/datasets/alignments/wes/v3/bam/ascat_trycatch
