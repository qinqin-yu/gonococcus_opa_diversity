This directory contains the Snakemake pipeline and downstream Jupyter notebooks to run the main analysis of *opa* diversity, phase, and evolution for the paper.

The Snakemake pipeline is compatible with Snakemake v8+ and conda v24.7.1+. 

## Installing snakemake
1. Make a new environment: `conda create --name snakemake8`
2. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html): `mamba create -c conda-forge -c bioconda -n snakemake snakemake=8.25.5` (can alternatively use conda)
3. [Install `executor-plugin-cluster-generic`](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html): `mamba install snakemake-executor-plugin-cluster-generic`

## To run
1. Deposit complete genomes into the directory `input_data/complete_genome_assemblies/`.
2. Submit job using `sbatch start_snakemake.sh` from the main directory.

## Manual scripts
The following parts of the analysis must be run manually:
- All Jupyter notebooks to analyze the outputs from the Snakemake pipeline (`workflow/notebooks/`)
- iTol tree must be generated in [iTol](https://itol.embl.de/) (`results/itol/complete_genome_pseudogenomes.final_tree.itol_order.tre`)
- BEAST analysis to get time-scaled phylogeny of subtree (Generate the file `results/beast/complete_genomes_recomb_masked.xml` in BEAUTI and then run `workflow/scripts/beast.sh`)
