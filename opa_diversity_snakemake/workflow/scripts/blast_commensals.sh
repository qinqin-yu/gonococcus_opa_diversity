#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4G
#SBATCH -p shared
#SBATCH -t 0-01:00
#SBATCH -o blast.out
#SBATCH -e blast.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --account=grad_lab

mkdir -p blastdb
mkdir -p blast
mkdir -p assemblies_tmp

# Replaces the header in each fasta file with the name of the file (excluding ".fa") and concatenates them into one fasta file

for file in assemblies/*.fa; do
    base=$(basename "$file" .fa)
    sed "1s/.*/>$base/" "$file" > "assemblies_tmp/$base.fa"
done

cat assemblies_tmp/*.fa > blastdb/gc_contigs.fa

rm -r assemblies_tmp/

# Make blast database

makeblastdb -dbtype nucl -in blastdb/gc_contigs.fa -out blastdb/gc

# Run blast
# -outfmt 7 gives a commented tab-delimited file
# max_hsps sets the maximum matches in a subject

blastn -db blastdb/gc -query FA1090_opa_1_reference_122_263.fa -out blast/opa_commensals_con1_blast.tab -num_threads 1 -max_hsps 15 -outfmt 7
blastn -db blastdb/gc -query FA1090_opa_1_reference_336_455.fa -out blast/opa_commensals_con2_blast.tab -num_threads 1 -max_hsps 15 -outfmt 7
blastn -db blastdb/gc -query FA1090_opa_1_reference_598_732.fa -out blast/opa_commensals_con3_blast.tab -num_threads 1 -max_hsps 15 -outfmt 7

