#!/bin/bash
#SBATCH -p shared
#SBATCH -c 12
#SBATCH -t 0-03:00
#SBATCH --mem-per-cpu=12G
#SBATCH -o gubbins_nanopore_genomes.out
#SBATCH -e gubbins_nanopore_genomes.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu 

run_gubbins.py --first-tree-builder rapidnj --tree-builder raxmlng --first-model JC --model GTR --threads 12 pseudogenome_alignment.fa
