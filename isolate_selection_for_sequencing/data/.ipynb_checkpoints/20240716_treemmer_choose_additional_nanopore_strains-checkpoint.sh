#!/bin/bash
#SBATCH -p shared
#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH --mem-per-cpu=12G
#SBATCH -o 20240716_treemmer_choose_additional_nanopore_strains.out
#SBATCH -e 20240716_treemmer_choose_additional_nanopore_strains.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu 

singularity exec /n/holylfs05/LABS/grad_lab/Lab/software/treemmer_sb/ python3 /Treemmer_v0.3.py --list_meta treemmer_metadata_input.csv --list_meta_count list_meta_count.txt --stop_at_X_leaves 210 --cpu 4 lab_strains_hybrid_genomes.final_tree.tre > RTL.txt