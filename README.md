# gonococcus_opa_diversity

This repository contains the code to reproduce the analysis in:

QinQin Yu, Tatum D. Mortimer, Sofia Blomqvist, Bailey Bowcutt, David Helekal, Samantha G. Palace, Yonatan H. Grad.  Diversity and evolution of a phase-variable multi-locus antigen in *Neisseria gonorrhoeae*.

The contents of the repository are described below: 

## complete_genome_assembly
This directory contains code to assemble complete genomes from long-read Nanopore sequencing data. It also contains code to test the effect of different assembly and polishing methods on the accuracy of *opa* sequences.

## draft_genome_assembly
This directory contains code to assemble draft genomes using short-read (Illumina) sequencing data. It also contains code to simulate short reads from complete genomes (to ensure that their draft genomes are assembled in the same way). 

## isolate_selection_for_sequencing
This directory contains the code, data, and figures for determining the isolates selected for Nanopore sequencing. The isolates were selected from across the *N. gonorrhoeae* phylogeny after accounting for the diversity that is represented in the publicly available complete genomes.

## opa_diversity_snakemake
This directory contains the Snakemake pipeline to run the main analysis for the paper assessing *opa* diversity, phase, and evolution. More details about the contents of the pipeline and how to run it are in the README file in that directory. 
