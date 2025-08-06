#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH -p shared
#SBATCH -t 0-01:00
#SBATCH -o mafft_mcl_clusters.out
#SBATCH -e mafft_mcl_clusters.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --account=grad_lab
#SBATCH --array=0-219

mkdir -p cluster_alignments/

readarray -t f < cluster_sequences/filenames.txt

mafft --auto cluster_sequences/${f[${SLURM_ARRAY_TASK_ID}]} > cluster_alignments/${f[${SLURM_ARRAY_TASK_ID}]}.aln