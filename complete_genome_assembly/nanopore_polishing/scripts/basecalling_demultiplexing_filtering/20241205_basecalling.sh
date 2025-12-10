#!/bin/bash
#SBATCH --job-name=basecalling              # Job name
#SBATCH --partition=gpu             # Partition (queue) name
#SBATCH --gres=gpu:1                  # Requests one GPU device 
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=2             # Number of CPU cores per task
#SBATCH --mem=20gb                    # Job memory request
#SBATCH --time=03:00:00               # Time limit hrs:min:sec
#SBATCH --output=basecalling.out         # Standard output log
#SBATCH --error=basecalling.err          # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu  # Where to send mail (change to your email address)

dorado basecaller sup /n/holyscratch01/grad_lab/Lab/nanopore_data_transfer/20240829_representative_genomes_run2/no_sample/20240829_1559_MN36583_FAZ94693_31d3c394/pod5 -r --kit-name SQK-NBD114-24 > calls.bam
