# ENA data submission

This directory contains the scripts and files for submitting the raw Nanopore sequencing reads and complete genome assemblies to the European Nucleotide Archive (ENA). 

## Raw reads submission

Steps:

1. Register a new study (project) in the ENA Webin portal. The study accession will be needed for submitting the reads.
2. For any isolates where there is not an existing biosample accession, register the samples in the ENA Webin portal. The biosample accessions will be needed for submitting the reads.
3. Collect the gzipped fastq read files into a folder. If the read files are in bam format, convert to gzipped fastq files using the snakemake pipeline in `reads_submission/scripts/convert_bam_to_fastq`.
4. Generate manifest (metadata) files for each submission (in this case, for each gzipped fastq file). Can create a file for all samples (see for example, `reads_submission/data/ena_reads_submission_manifest_combined.tsv`) and then separate into individual files using the script in `reads_submissions/scripts/generate_individual_manifest_files.ipynb`). The individual manifest files can be found in `reads_submission/data/reads_manifest_files/`. More info about ENA manifest files can be found [here](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html) and [here](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html#webin-cli-submission).
5. Validate and then submit the reads and manifest files using the snakemake pipeline in `reads_submission/scripts/ena_web_cli_read_submission`. You will need to get the Grad lab ENA Webin account login information from Grad lab manager.

Notes:
- The most raw sequencing reads (basecalled and demultiplexed) were submitted.
- In this case, almost all biosamples (the data on the isolates themselves) had been previously submitted to ENA or NCBI when they were originally published (when they were sequenced with short read sequencing). In this case, the recommended best practice when uploading new data (in this case, long read sequencing data) is to link the new data with the existing biosample accessions. NCBI does not let you do this if the biosample was originally submitted to ENA, which is what motivated me to submit using ENA for this project. Additionally, if the biosample was submitted many years ago, it's possible that some of the fields have invalid values, preventing the linking to the new data. In this case, my workaround was just to create a new biosample for these particular isolates.

## Assembly submission
The process for submitting the complete genomes is similar as for submitting the reads, just with an additional "chromosome list" file.

Steps:

1. Collect the gzipped fasta read files into a folder. 
2. Generate manifest (metadata) files for each submission (in this case, for each gzipped fasta file). Can create a file for all samples (see for example, `assembly_submission/data/ena_assembly_submission_manifest_combined.tsv`) and then separate into individual files using the script in `assembly_submissions/scripts/generate_individual_manifest_files.ipynb`). The individual manifest files can be found in `assembly_submission/data/assembly_manifest_files/`. More info about ENA manifest files can be found [here](https://ena-docs.readthedocs.io/en/latest/submit/assembly/genome.html) and [here](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html#webin-cli-submission).
3. Generate chromosome list files for each subimssion using `assembly_submission/scripts/write_chromosome_list_files.ipynb`. More info about ENA chromosome list files can be found [here](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/assembly.html#chromosome-list-file).
4. Validate and then submit the reads and manifest files using the snakemake pipeline in `assembly_submission/scripts/ena_web_cli_assembly_submission`. You will need to get the Grad lab ENA Webin account login information from Grad lab manager.

## Making bioproject public

The bioproject can be made public on the ENA online Webin portal upon manuscript submission. Note that it took about 2 days for everything to appear, so I did this about 2 days before submission. 