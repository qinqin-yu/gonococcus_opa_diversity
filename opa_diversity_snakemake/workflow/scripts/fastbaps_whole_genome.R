#!/usr/bin/env Rscript

library(fastbaps)
library(ggtree)
library(phytools)
library(ggplot2)
library(argparse)

parser <- ArgumentParser(description='Run fastbaps to cluster whole genomes')
parser$add_argument('--pseudogenome_alignment', help='Whole genome pseudogenome alignment')
parser$add_argument('--tree', help='Gubbins whole genome tree')
parser$add_argument('--clusters', help='Output clusters')

args<-parser$parse_args()

# Partition based on tree
sparse.data <- import_fasta_sparse_nt(args$pseudogenome_alignment, prior = "baps")

# Note that the final_tree.tre output of gubbins contains :0.0; at the end of the file. Need to delete :0.0 (keeping the semicolon) for the tree loading to work. 

gubbins.tree <- phytools::read.newick(args$tree)
gubbins.rooted <- phytools::midpoint.root(gubbins.tree)
best.partition <- best_baps_partition(sparse.data, gubbins.rooted)
best.partition.df <- data.frame(id = gubbins.rooted$tip.label, fastbaps = best.partition, stringsAsFactors= FALSE)
write.csv(best.partition.df, paste(args$clusters, sep = ""), row.names=FALSE)