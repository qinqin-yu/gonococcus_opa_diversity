#!/bin/bash
#SBATCH --job-name=filter              # Job name
#SBATCH --partition=hsph             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                    # Job memory request
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --output=filter.out         # Standard output log
#SBATCH --error=filter.err          # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu  # Where to send mail (change username@uga.edu to your email address)
#SBATCH --array=0-15

readarray -t f < barcode_numbers.txt

mkdir -p filtered_fastqs

filtlong -p 90 --min_length 1000 demultiplexed_reads/SQK-NBD114-24_barcode${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz | gzip > filtered_fastqs/barcode${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz
