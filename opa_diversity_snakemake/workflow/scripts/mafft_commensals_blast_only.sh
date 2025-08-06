#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH -p shared
#SBATCH -t 0-01:00
#SBATCH -o mafft.out
#SBATCH -e mafft.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --account=grad_lab

mafft --auto --op 4 --ep 1 opa_blast_con1_only_with_MC58.fa > opa_blast_con1_only_with_MC58.fa.aln
mafft --auto --op 4 --ep 1 opa_blast_con3_only_with_MC58.fa > opa_blast_con3_only_with_MC58.fa.aln