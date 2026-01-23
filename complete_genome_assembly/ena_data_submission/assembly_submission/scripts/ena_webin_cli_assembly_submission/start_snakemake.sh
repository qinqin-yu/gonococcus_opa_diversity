#!/bin/bash
#SBATCH --job-name=simulate_short_reads
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH	-p shared
#SBATCH -t 1-00:00
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH --account=grad_lab

mkdir -p outputs

snakemake -k --profile config --latency-wait 10 --cores 12
