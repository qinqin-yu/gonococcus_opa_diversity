#!/bin/bash
#SBATCH --job-name=opa_diversity
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH	-p shared
#SBATCH -t 1-00:00
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH --account=grad_lab

snakemake -k --profile config --latency-wait 10 --omit-from gubbins
