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

mafft --amino excluded_opa_genes_aa.fa > excluded_opa_genes_aa.fa.aln
mafft --nuc excluded_opa_genes_nt.fa > excluded_opa_genes_nt.fa.aln
mafft --nuc excluded_opa_genes_nt_type1.fa > excluded_opa_genes_nt_type1.fa.aln
mafft --nuc excluded_opa_genes_nt_type2.fa > excluded_opa_genes_nt_type2.fa.aln
mafft --nuc excluded_opa_genes_nt_type3.fa > excluded_opa_genes_nt_type3.fa.aln
mafft --nuc excluded_opa_genes_nt_type4.fa > excluded_opa_genes_nt_type4.fa.aln
