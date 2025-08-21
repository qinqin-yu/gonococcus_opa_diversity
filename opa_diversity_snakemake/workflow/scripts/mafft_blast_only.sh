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

mafft --nuc opa_con1_blast_only_with_FA1090.fa > opa_con1_blast_only_with_FA1090.fa.aln
mafft --nuc opa_con2_blast_only_with_FA1090.fa > opa_con2_blast_only_with_FA1090.fa.aln
mafft --nuc opa_con3_blast_only_with_FA1090.fa > opa_con3_blast_only_with_FA1090.fa.aln