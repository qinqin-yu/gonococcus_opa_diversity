#!/bin/bash
#SBATCH -p shared
#SBATCH -c 12
#SBATCH -t 0-06:00
#SBATCH --mem-per-cpu=12G
#SBATCH -o 20250219_beast_samples_with_dates.out
#SBATCH -e 20250219_beast_samples_with_dates.err
#SBATCH --mail-type=END
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu 
#SBATCH --array=1-7

/n/holylfs05/LABS/grad_lab/Users/qinqinyu/software/beast/bin/beast -threads 12 subtree_${SLURM_ARRAY_TASK_ID}_recomb_masked_with_dates.xml
