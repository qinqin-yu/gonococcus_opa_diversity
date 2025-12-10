## autocycler_assembly
Contains the Snakemake pipeline to reproduce the assembly of complete genomes from the manuscript using Nanopore long-read sequencing data. This pipeline uses the program [Autocycler](https://github.com/rrwick/Autocycler). More information can be found on the README file within the directory.

## autocycler_opa_accuracy
Contains the code used to test the effect of sequencing read depth on the *opa* sequences from the Autocycler assemblies.

## nanopore_polishing
Contains the code used to test the effect of different assembly and polishing methods on the *opa* sequences. The polishing methods tested include those using short-read and long-read sequencing data.
