#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio import Phylo

import pairwise_distance_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Calculate between-strain opa repertoire distance')
    parser.add_argument("alignment_filename", help="Opa amino acid alignment filename")
    parser.add_argument("phylogeny", help="itol complete genomes tree in newick format")
    parser.add_argument("pairsnp_output", help="Pairsnp output of pseudogenome distances")
    parser.add_argument("distance_df", help="csv of distances by matched opa")
    parser.add_argument("opa_genome_distance", help="csv of opa and genomic distance between strains")
    return parser.parse_args()

args = get_args()

# Calculate the distance between opa repertoires in the amino acid alignment
filename = args.alignment_filename
aln = AlignIO.read(filename, 'fasta')
distance_df, summary = fnc.opa_repertoire_distance(aln)
distance_df.to_csv(args.distance_df)

# Calculate the distances between genomes and add to summary df

# Load the whole genome tree
trees = Phylo.parse(args.phylogeny, "newick")
for tree in trees:
    break
    
# Calculate the whole genome distance between strains (branch length on the tree) and save to df

genome_distances = []
for i, row in summary.iterrows():
    strain_A = row['strain_A']
    strain_B = row['strain_B']
    genome_distances.append(tree.distance(strain_A.replace('#', '_'), strain_B.replace('#', '_')))
summary['genome_distance']=genome_distances

# Also add pseudogenome distances

genome_distance = pd.read_csv(args.pairsnp_output, header = None)

names = genome_distance[0]
distances = genome_distance.drop(0, axis = 'columns').values

for i, row in summary.iterrows():
    strain_A = row['strain_A']
    strain_B = row['strain_B']
    idx_A = np.where(names==strain_A)[0][0]
    idx_B = np.where(names==strain_B)[0][0]
    summary.at[i, 'pseudogenome_distance'] = distances[idx_A, idx_B]

summary.to_csv(args.opa_genome_distance)