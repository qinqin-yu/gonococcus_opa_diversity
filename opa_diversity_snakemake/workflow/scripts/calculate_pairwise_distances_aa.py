#!/usr/bin/env python

import glob
import os
import argparse
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

from scipy.stats import kstest

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
font = {'family' : 'Arial',
        'size'   : 12}
matplotlib.rc('font', **font)

import pairwise_distance_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Calculate pairwise distances and plot results')
    parser.add_argument("alignment_filename", help="Opa amino acid alignment filename")
    parser.add_argument("pairwise_distance_filename", help="Pairwise distance output csv filename")
    parser.add_argument("png_filename", help="png filename of figure comparing within and between strain opa diversity")
    parser.add_argument("pdf_filename", help="pdf filename of figure comparing within and between strain opa diversity")
    return parser.parse_args()

### CALCULATE PAIRWISE DISTANCES 

# Load nucleotide alignment

args = get_args()
filename =  args.alignment_filename
aln = AlignIO.read(filename, 'fasta')

# Calculate distance between all pairs of sequences

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

names = np.array(dm.names)
matrix = dm.matrix

# Convert distance output list to array

matrix_array = np.empty((len(matrix), len(matrix)))
matrix_array.fill(np.nan)
for i in range(len(matrix)):
    matrix_array[i, 0:i+1] = matrix[i]

# Replace zeros on diagonal with nan
np.fill_diagonal(matrix_array, np.nan)

matrix_array_flattened = matrix_array.flatten()
matrix_array_flattened_nonan = matrix_array_flattened[~np.isnan(matrix_array_flattened)]

# Get a dataframe of the within-strain opa pairwise distances

# First get the strain name from the opa gene name

ids = []
for i in range(len(aln)):
    ids.append(aln[i].id)
    
strains = []
opa_gene_nums = []
for id_i in ids:
    strain, opa_gene_num = fnc.parse_opa_name(id_i)
    strains.append(strain)
    opa_gene_nums.append(opa_gene_num)

ids_metadata = pd.DataFrame({'id':ids, 'strain':strains, 'opa_gene_nums':opa_gene_nums})

# Then get the pairwise distances between opas within each strain

within_strain_distance = pd.DataFrame()
for strain, df in ids_metadata.groupby('strain'):
    idx_strain = df.index
    aln_strain = aln[np.min(idx_strain):np.max(idx_strain)+1]
    matrix_strain_array, distance_df = fnc.get_distance_matrix(aln_strain, get_df = True)
    within_strain_distance = pd.concat([within_strain_distance, distance_df])
within_strain_distance.reset_index(inplace = True, drop = True)

# Save df
within_strain_distance.to_csv(args.pairwise_distance_filename)

# Calculate whether within-strain opa distances are significantly different than that between all strains (within and between strains)
strains = []
pvalues = []
for strain, df in within_strain_distance.groupby('strain_A'):
    result = kstest(df['distance'], matrix_array_flattened_nonan)
    strains.append(strain)
    pvalues.append(result.pvalue)
    
compare_within_between_dist = pd.DataFrame({'strain':strains, 'p_value':pvalues})
compare_within_between_dist['q_value']=compare_within_between_dist['p_value']*len(compare_within_between_dist)
sig_strains = compare_within_between_dist[compare_within_between_dist['q_value']<0.05]['strain'].values

# Get a dataframe of the between strain opa pairwise distances

# First get the strain name from the opa gene name

matrix_array_between_strains = np.copy(matrix_array)
# between_strain_distance = pd.DataFrame()
for strain, df in ids_metadata.groupby('strain'):
    idx_strain = df.index
    matrix_array_between_strains[np.min(idx_strain):np.max(idx_strain)+1, np.min(idx_strain):np.max(idx_strain)+1] = np.nan

between_strains_flattened = matrix_array_between_strains.flatten()
between_strains_flattened_nonan = between_strains_flattened[~np.isnan(between_strains_flattened)]

### PLOT RESULTS

# Plot the comparison of pairwise distances within strains and between strains

# plt.figure(figsize = (6,3))
# bins = np.linspace(0, 0.5, 50)
# plt.hist(between_strains_flattened_nonan, bins = bins, label = 'Between strains', density = True, histtype = 'step')
# plt.hist(within_strain_distance['distance'].values, bins = bins, label = 'Within strains', density = True, histtype = 'step')
# plt.xlabel('Pairwise amino acid distance across $opa$ genes')
# plt.ylabel('Density')
# plt.legend(loc = (1.01, 0))
# plt.yscale('log')
# plt.tight_layout()
# plt.savefig(args.png_filename, dpi = 300)
# plt.savefig(args.pdf_filename)
# plt.show()

plt.figure(figsize = (6,3))
bins = np.linspace(0, 0.5, 50)
plt.hist(between_strains_flattened_nonan, bins = bins, label = 'Between isolates', density = True, histtype = 'step', color = '#919BCA')
plt.hist(within_strain_distance['distance'].values, bins = bins, label = 'Within isolates', density = True, histtype = 'step', color = '#ad6a6c')
plt.xlabel('Pairwise amino acid distance\nbetween $opa$ genes')
plt.ylabel('Density')
plt.legend(loc = (1.01, 0))
plt.yscale('log')
plt.tight_layout()
plt.savefig(args.png_filename, dpi = 300)
plt.savefig(args.pdf_filename)
plt.show()