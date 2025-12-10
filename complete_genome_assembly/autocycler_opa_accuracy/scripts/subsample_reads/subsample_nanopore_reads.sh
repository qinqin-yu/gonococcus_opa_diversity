#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH	-p shared
#SBATCH -t 0-01:00
#SBATCH -o subsample_nanopore_reads.out
#SBATCH -e subsample_nanopore_reads.err
#SBATCH --account=grad_lab
#SBATCH --array=0-5

readarray -t f < samples.txt

mkdir -p subsampled_reads/

seqtk sample -s100 reads/${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz 0.15 > subsampled_reads/${f[${SLURM_ARRAY_TASK_ID}]}_15.fastq
seqtk sample -s100 reads/${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz 0.3 > subsampled_reads/${f[${SLURM_ARRAY_TASK_ID}]}_30.fastq
seqtk sample -s100 reads/${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz 0.6 > subsampled_reads/${f[${SLURM_ARRAY_TASK_ID}]}_60.fastq
gzip subsampled_reads/*.fastq