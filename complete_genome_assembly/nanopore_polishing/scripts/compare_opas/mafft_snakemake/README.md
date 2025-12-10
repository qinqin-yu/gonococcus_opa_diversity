# Pipeline for annotating opa genes in complete *N. gonorrhoeae* and *N. meningitidis* genomes

## Pipeline description
This pipeline identifies the location of opa genes in complete genomes and saves their sequences. 

## Installing snakemake
1. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) into a new environment

## Running the pipeline
1. Make a directory for your current project
2. Copy the contents of this directory (`Snakefile`, `start_snakemake.sh`, `scripts/`, `conda_envs/`, `slurm_genomics_pipeline/`)
3. Edit the snakemake config (`slurm_genomics_pipeline/config.yaml`) with your conda installation path, specifically the line that starts with `conda-prefix:`
4. Edit the submit script (`start_snakemake.sh`) with your email address (you can also update the time/queue/etc. depending on how many genomes you have)
5. Add executable permissions to `slurm-status.py` and `scripts/*.py` if needed (`chmod +x slurm_genomics_pipeline/slurm-status.py`, `chmod +x scripts/*.py`)
6. Make a directory called `assemblies` with all your complete genome assemblies with the format `[sample_name].fa`. If you have used the [Nanopore_analysis](https://github.com/gradlab/Nanopore_analysis) pipeline, you can use this [script](https://hu.sharepoint.com/sites/GradLab/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2FGradLab%2FShared%20Documents%2FBailey%2Frename%5Ffasta%5Fmove%5Fto%5Fopa%5Fannotations%2Esh&parent=%2Fsites%2FGradLab%2FShared%20Documents%2FBailey) to rename the consensus sequence files. Make sure there are no spaces in the filenames. Also the assemblies must be complete (cannot have multiple contigs).
7. Activate the environment that has snakemake installed.
8. Submit job using `sbatch start_snakemake.sh`

## Pipeline outputs

```
opa_sequences/
    A fasta file for each sample containing the sequence in the open reading frame of each opa gene
opa_locations/
    A csv file for each sample containing the location in the genome of the open reading frame, coding repeat sequences, termination sequence, and strand on which each opa gene is found
```
## Notes
`opa_locations/*.csv` also keeps track of truncated genes where only the sequence near the stop codon has been found. These truncated genes are not given an identifier (the `id` field is blank) and you can use this to remove these from downstream analyses. 
