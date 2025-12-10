# Nanopore_analysis
Scripts and instructions for Nanopore data processing and analysis

## Assembly pipeline
Takes **BASECALLED NANOPORE READS** as input, assembles with [Autocycler](https://github.com/rrwick/Autocycler/wiki). Autocycler is an ensemble method that subsets the read dataset, assemblies the read subsets using different assembly methods, and combines the different assemblies to create a consensus assembly. It's been shown that Autocycler generally produces more accurate bacterial assemblies than any assembler on its own. This pipeline requires read depths of at least 25x and recommended read depths of 50x.

This pipeline is for **Nanopore-only** assemblies. In most cases, Nanopore-only assemblies created with Autocycler will be sufficient, but in some cases when you want to remove any remaining errors you may want to consider additional short-read polishing (see notes below).

The code is adapted from the [Autocycler wiki](https://github.com/rrwick/Autocycler/wiki/Fully-automated-assembly). 

### Required input

In working directory:

* Subdirectory named `fastq/` containing basecalled, demultiplexed reads (usually from dorado). The fastq files should be gzipped so that they have file names {sample_name}.fastq.gz
* `slurm_genomics_pipeline/` subdirectory: contains config information for snakemake
* `conda_envs/` subdirectory: contains YAML files for making Autocycler conda environment which includes the dependencies of the individual assembly programs
* `scripts/` subdirectory: Autocycler program and helper scripts for running individual assembly programs
* `Snakefile` containing pipeline instructions
* `start_snakemake.sh` script to submit the job to the cluster. Note that this file must be executable! To change permissions, use the command `chmod +x start_snakemake.sh`

### Expected output

In working directory:

* `fastqs/` subdirectory: for each sample, will contain a single file with all input fastq reads concatenated into a single file and gzipped
* `autocycler/` subdirectory: all output
	* `subsampled_reads/` subdirectory: 4 sets of subsampled reads. 
	* `assemblies/` subdirectory: assemblies created using different programs for each of the subsampled read sets
	* `autocycler_out/` subdirectory: outputs of running Autocycler to combine the assemblies. The final consensus assemblies are at `autocycler_out/consensus_assembly.fasta`
* `consensus_assmeblies/` The final consensus assemblies copied from `autocycler_out/{sample}/consensus_assembly.fasta`
* `logs/slurm/` subdirectory: snakemake log files
* `metrics.tsv` file: tab-delimited file containing metrics from the autocycler assembly to assess quality. The description of the metrics can be found on the [Autocycler wiki](https://github.com/rrwick/Autocycler/wiki/Autocycler-table#default-fields).
* `autocycler.err` file: contains everything that the processes in the Snakefile pipe to STDOUT during the pipeline's run. You can check this as the pipeline runs to see where you are. You can also look at it if something fails to try to figure out what went wrong.
* `autocycler.out` file: contains anything piped to STDOUT when the cluster submission finishes; usually this file is empty if all has gone well!

In the directory where your conda environments are stored:

The pipeline will install a new Autocycler environment for its own use. These will be named some crazy long string of characters and will not conflict with existing Autocycler environments you might have. You can manually delete these once you are satisfied that you have the output you are looking for.

### Running the pipeline

#### Edit the following files:

* In `slurm_genomics_pipeline/config.yaml`, change the path listed next to `conda-prefix:` to the directory where your own conda environments are stored. (This is where the flye and medaka environments for the pipeline will be installed.)
* In `start_snakemake.sh`, enter your own email address in the header to receive an email when the job finishes successfully. Optionally, you can add additional #SBATCH instructions to e.g. get emails if the job quits as a result of an error.

#### Submit the job

Activate the snakemake conda environment. (This does not come pre-installed on the cluster, so you may need to install it first.) Note that this pipeline uses snakemake versions that are 7 or earlier (will not work with snakemake 8).

Submit your job to the cluster using the following command: `sbatch start_snakemake.sh`

You can check on the progress of your job with `squeue`, `sacct`, or by looking at the contents of `autocycler.err`.

### Quality control
It's a good idea to check the quality of the consensus assemblies after the run has finished. A summary of metrics of the input data and assemblies are stored in the `metrics.tsv` file. For basic quality control, you can check that the consensus assembly was fully resolved (i.e. one complete genome plus optional plasmids) in the column named `consensus_assembly_fully_resolved` and that the consensus assembly length is about what you would expect for GC (~2.2 Mbp) in the column named `consensus_assembly_bases`. More information on the other metrics of the table and other checks for quality control can be found on the [Autocycler wiki](https://github.com/rrwick/Autocycler/wiki/Autocycler-table#default-fields). 

### Notes

* The currently recommendations are to not polish with medaka. In most instances, Autocycler-only assemblies will be sufficient, but polishing with Polypolish default and Pypolca careful can potentially fix any remaining errors in the assembly (1). 
* This pipeline only works with read depths of at least 25x. The recommendation is to aim for read depths of 50x. If you have lower read depth than 25x, it is possible to modify the pipeline to run with lower read depth, but it's not clear what the accuracy of the assemblies will be.
  
