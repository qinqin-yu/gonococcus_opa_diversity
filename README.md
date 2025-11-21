# gonococcus_opa_diversity

This repository contains the code to reproduce the analysis in:

QinQin Yu, Tatum D. Mortimer, Sofia Blomqvist, Bailey Bowcutt, David Helekal, Samantha G. Palace, Yonatan H. Grad.  Diversity and evolution of a phase-variable multi-locus antigen in *Neisseria gonorrhoeae*.

The contents of the repository are described below: 

## draft_genome_assembly
This directory contains code to assemble draft genomes using short-read (Illumina) sequencing data. It also contains code to simulate short reads from complete genomes (to ensure that their draft genomes are assembled in the same way). 

## opa_diversity_snakemake
This directory contains the Snakemake pipeline to run the main analysis for the paper assessing *opa* diversity, phase, and evolution. More details about the contents of the pipeline and how to run it are in the README file in that directory. 
