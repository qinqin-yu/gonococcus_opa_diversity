#!/bin/bash
#SBATCH --job-name=autocycler              # Job name
#SBATCH --partition=shared             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task	
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=20gb                    # Job memory request
#SBATCH --time=06:00:00               # Time limit hrs:min:sec
#SBATCH --output=autocycler.out         # Standard output log
#SBATCH --error=autocycler.err          # Standard error log
#SBATCH --array=0-17

#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=qinqinyu@hsph.harvard.edu  # Where to send mail (change to your email address)

readarray -t f < samples.txt

mkdir -p autocycler

reads=fastqs/${f[${SLURM_ARRAY_TASK_ID}]}.fastq.gz  # your read set goes here
threads=16  # set as appropriate for your system
genome_size=$(genome_size_raven.sh "$reads" "$threads")  # can set this manually if you know the value

# Step 1: subsample the long-read set into multiple files
autocycler subsample --reads "$reads" --out_dir autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/subsampled_reads --genome_size "$genome_size"

# Step 2: assemble each subsampled file
mkdir autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/assemblies
for assembler in canu flye miniasm necat nextdenovo raven; do
    for i in 01 02 03 04; do
        "$assembler".sh autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/subsampled_reads/sample_"$i".fastq autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/assemblies/"$assembler"_"$i" "$threads" "$genome_size"
    done
done

# Optional step: remove the subsampled reads to save space
rm autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/assemblies -a autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/autocycler_out

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/autocycler_out

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    autocycler resolve -c "$c"
done

# Step 7: combine resolved clusters into a final assembly
autocycler combine -a autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/autocycler_out -i autocycler/${f[${SLURM_ARRAY_TASK_ID}]}/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa
