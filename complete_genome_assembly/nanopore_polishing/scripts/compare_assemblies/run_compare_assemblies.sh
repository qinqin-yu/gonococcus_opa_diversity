#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH-p shared
#SBATCH -t 0-00:10
#SBATCH -o compare_assemblies.out
#SBATCH -e compare_assemblies.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu
#SBATCH --array=0-5

readarray -t f < accessions.txt

flye_path="../../data/kit14/flye_assemblies/"
medaka_path="../../data/kit14/medaka_polished_assemblies/"
polypolish_path="../../data/kit14/polypolish/"
pypolca_path="../../data/kit14/pypolca/"

flye_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241120_compare_assemblies/flye_assemblies_sorted/"
medaka_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241120_compare_assemblies/medaka_polished_assemblies_sorted/"
polypolish_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241120_compare_assemblies/polypolish_sorted/"
pypolca_sorted_path="/n/netscratch/grad_lab/Lab/qinqinyu/20241120_compare_assemblies/pypolca_sorted/"

mkdir -p $flye_sorted_path
mkdir -p $medaka_sorted_path
mkdir -p $polypolish_sorted_path
mkdir -p $pypolca_sorted_path

# First sort contigs so that they are in the same order in each fasta

seqkit sort -n $flye_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta>$flye_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta
seqkit sort -n $medaka_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta>$medaka_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta
seqkit sort -n $polypolish_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.fasta>$polypolish_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.fasta
seqkit sort -n $pypolca_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.pypolca.fasta>$pypolca_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.pypolca.fasta

# Compare flye and medaka

flye_medaka_out="../../data/kit14/compare_assemblies/flye_medaka/"
mkdir -p $flye_medaka_out

./compare_assemblies.py $flye_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta $medaka_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta>$flye_medaka_out${f[${SLURM_ARRAY_TASK_ID}]}.txt

grep -o "*" $flye_medaka_out${f[${SLURM_ARRAY_TASK_ID}]}.txt | wc -l > $flye_medaka_out${f[${SLURM_ARRAY_TASK_ID}]}.num_diff.txt

# Compare medaka and polypolish

medaka_polypolish_out="../../data/kit14/compare_assemblies/medaka_polypolish/"
mkdir -p $medaka_polypolish_out

./compare_assemblies.py $medaka_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.fasta $polypolish_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.fasta>$medaka_polypolish_out${f[${SLURM_ARRAY_TASK_ID}]}.txt

grep -o "*" $medaka_polypolish_out${f[${SLURM_ARRAY_TASK_ID}]}.txt | wc -l > $medaka_polypolish_out${f[${SLURM_ARRAY_TASK_ID}]}.num_diff.txt

# Compare polypolish and pypolca

polypolish_pypolca_out="../../data/kit14/compare_assemblies/polypolish_pypolca/"
mkdir -p $polypolish_pypolca_out

./compare_assemblies.py $polypolish_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.fasta $pypolca_sorted_path${f[${SLURM_ARRAY_TASK_ID}]}.polypolish.pypolca.fasta>$polypolish_pypolca_out${f[${SLURM_ARRAY_TASK_ID}]}.txt

grep -o "*" $polypolish_pypolca_out${f[${SLURM_ARRAY_TASK_ID}]}.txt | wc -l > $polypolish_pypolca_out${f[${SLURM_ARRAY_TASK_ID}]}.num_diff.txt