#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH	-p shared
#SBATCH -t 1-00:00
#SBATCH -o nanopore_flye_medaka.out
#SBATCH -e nanopore_flye_medaka.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu
#SBATCH --account=grad_lab

mkdir -p logs/slurm
snakemake -k --rerun-incomplete --profile slurm_genomics_pipeline --set-resources medaka:slurm_partition=gpu
