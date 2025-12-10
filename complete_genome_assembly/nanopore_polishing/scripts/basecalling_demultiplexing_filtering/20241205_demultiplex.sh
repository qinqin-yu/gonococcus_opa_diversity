#!/bin/bash
#SBATCH --job-name=demux              # Job name
#SBATCH --partition=gpu             # Partition (queue) name
#SBATCH --gres=gpu:1                  # Requests one GPU device 
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=2             # Number of CPU cores per task
#SBATCH --mem=10gb                    # Job memory request
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH --output=demux.out         # Standard output log
#SBATCH --error=demux.err          # Standard error log

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu  # Where to send mail (change to your email address)

dorado demux --output-dir demultiplexed_reads --no-classify --emit-fastq calls.bam
gzip demultiplexed_reads/*.fastq
rename "9095f4d123708578d426730579a5a748b9ba3fb5_" "" demultiplexed_reads/*.fastq.gz
