#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH-p shared
#SBATCH -t 0-00:10
#SBATCH -o compare_assemblies.out
#SBATCH -e compare_assemblies.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu

pypolca_path="../../data/kit14/pypolca/"
pypolca_parent_with_evolved_path="../../data/kit14/pypolca_parent_with_evolved/"

pypolca_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241121_compare_assemblies/pypolca_sorted/"
pypolca_parent_with_evolved_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241121_compare_assemblies/pypolca_parent_with_evolved_sorted/"

mkdir -p $pypolca_sorted_path
mkdir -p $pypolca_parent_with_evolved_sorted_path

# First sort contigs so that they are in the same order in each fasta

seqkit sort -n ${pypolca_path}DGI_65.polypolish.pypolca.fasta>${pypolca_sorted_path}DGI_65.polypolish.pypolca.fasta
seqkit sort -n ${pypolca_path}EEE029.polypolish.pypolca.fasta>${pypolca_sorted_path}EEE029.polypolish.pypolca.fasta

seqkit sort -n ${pypolca_parent_with_evolved_path}DGI_65_a.polypolish.pypolca.fasta>${pypolca_parent_with_evolved_sorted_path}DGI_65_a.polypolish.pypolca.fasta
seqkit sort -n ${pypolca_parent_with_evolved_path}DGI_65_b.polypolish.pypolca.fasta>${pypolca_parent_with_evolved_sorted_path}DGI_65_b.polypolish.pypolca.fasta
seqkit sort -n ${pypolca_parent_with_evolved_path}EEE029_a.polypolish.pypolca.fasta>${pypolca_parent_with_evolved_sorted_path}EEE029_a.polypolish.pypolca.fasta
seqkit sort -n ${pypolca_parent_with_evolved_path}EEE029_b.polypolish.pypolca.fasta>${pypolca_parent_with_evolved_sorted_path}EEE029_b.polypolish.pypolca.fasta

# Compare polishing of parent with evolved

parent_with_evolved_out="../../data/kit14/compare_assemblies/parent_with_evolved/"
mkdir -p $parent_with_evolved_out

./compare_assemblies.py ${pypolca_sorted_path}DGI_65.polypolish.pypolca.fasta ${pypolca_parent_with_evolved_sorted_path}DGI_65_a.polypolish.pypolca.fasta>${parent_with_evolved_out}DGI_65_a.txt
./compare_assemblies.py ${pypolca_sorted_path}DGI_65.polypolish.pypolca.fasta ${pypolca_parent_with_evolved_sorted_path}DGI_65_b.polypolish.pypolca.fasta>${parent_with_evolved_out}DGI_65_b.txt
./compare_assemblies.py ${pypolca_sorted_path}EEE029.polypolish.pypolca.fasta ${pypolca_parent_with_evolved_sorted_path}EEE029_a.polypolish.pypolca.fasta>${parent_with_evolved_out}EEE029_a.txt
./compare_assemblies.py ${pypolca_sorted_path}EEE029.polypolish.pypolca.fasta ${pypolca_parent_with_evolved_sorted_path}EEE029_b.polypolish.pypolca.fasta>${parent_with_evolved_out}EEE029_b.txt

grep -o '*' ${parent_with_evolved_out}DGI_65_a.txt | wc -l > ${parent_with_evolved_out}DGI_65_a.num_diff.txt
grep -o '*' ${parent_with_evolved_out}DGI_65_b.txt | wc -l > ${parent_with_evolved_out}DGI_65_b.num_diff.txt
grep -o '*' ${parent_with_evolved_out}EEE029_a.txt | wc -l > ${parent_with_evolved_out}EEE029_a.num_diff.txt
grep -o '*' ${parent_with_evolved_out}EEE029_b.txt | wc -l > ${parent_with_evolved_out}EEE029_b.num_diff.txt
