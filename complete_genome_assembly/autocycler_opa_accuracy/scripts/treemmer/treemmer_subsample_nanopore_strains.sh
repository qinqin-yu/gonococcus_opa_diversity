#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH --mem-per-cpu=12G
#SBATCH -o treemmer_subsample_nanopore_strains.out
#SBATCH -e treemmer_subsample_nanopore_strains.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu 

singularity exec /n/holylfs05/LABS/grad_lab/Lab/software/treemmer_sb/ python3 /Treemmer_v0.3.py --list_meta treemmer_metadata_input.csv --list_meta_count list_meta_count.txt --stop_at_X_leaves 10 --cpu 4 nanopore_genomes.final_tree.tre > RTL.txt