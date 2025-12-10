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
#SBATCH --array=0-23

readarray -t f < barcode_numbers.txt

# Changes for future iterations of script:
# 1. Rename the fastqs as the wgs id of the strain - will need as input a mapping between the barcode numbers and the wgs ids.

samtools bam2fq -c 6 20240829_representative_genomes_run2_reads/9095f4d123708578d426730579a5a748b9ba3fb5_SQK-NBD114-24_barcode${f[${SLURM_ARRAY_TASK_ID}]}.bam | gzip > fastqs/barcode${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz

filtlong -p 90 --min_length 1000 fastqs/barcode${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz > filtered_fastqs/barcode${f[${SLURM_ARRAY_TASK_ID}]}_filtered.fastq.gz
