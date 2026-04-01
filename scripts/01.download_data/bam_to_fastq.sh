#!/bin/bash
#SBATCH -J gbm_download
#SBATCH -A FRANKELL-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --time=12:00:00
#SBATCH --array=1-75
#SBATCH -o logs/%A_%a_bam_to_fq.log
#SBATCH --mail-type=ALL

module load samtools

# get srr id
srr=$( cat srr.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
bam=${srr}.bam

echo "Processing ${bam}"

# index
samtools index bams/${bam}

# sort by name
echo "Sorting..."
samtools sort -@ 8 -n -o bams/${srr}.tmp.bam bams/${bam}

# convert to fastq
echo "Converting to fastqs..."
samtools fastq \
  -F 0x900 \
  -1 ${srr}_R1.fastq \
  -2 ${srr}_R2.fastq \
  -0 /dev/null \
  -s /dev/null \
  -n \
  bams/${srr}.tmp.bam

# zip
echo "Gzipping..."
gzip ${srr}_R1.fastq
gzip ${srr}_R2.fastq

# remove temporary bam
echo "Cleanup..."
rm bams/${srr}.tmp.bam

echo "Complete!"
