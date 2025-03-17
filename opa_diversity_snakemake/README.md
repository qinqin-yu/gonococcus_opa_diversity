Snakemake pipeline to run opa diversity and evolution analysis. Some parts of the analysis will be separate due to convenience (still working on figuring out which parts).

Compatible with Snakemake v8+ and conda v24.7.1+. You can install conda into the snakemake environment without modifying your main conda installation.

## Installing snakemake
1. Make a new environment (i.e. snakemake8)
2. [Install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
3. [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
4. [Install `executor-plugin-cluster-generic`](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html)

## To run
1. Deposit complete genomes into `input_data/complete_genome_assemblies/`.
2. Submit job using `sbatch start_snakemake.sh` from the main directory.

## Manual scripts
The following parts of the analysis must be run manually:
- All Jupyter notebooks (`workflow/notebooks/`)
- iTol tree (`results/itol/complete_genome_pseudogenomes.final_tree.itol_order.tre`)
- Beauti and BEAST2 (`results/beast/complete_genomes_recomb_masked.xml` and `workflow/scripts/beast.sh`)

## Notes
I've used the autocycler assembly for the complete genomes, as I've found that this is the most accurate assembly method.
