#!/usr/bin/env python
import argparse
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description='Assign opa gene loci')
    parser.add_argument("metadata_in", help="opa metadata from progressive mauve before assigning loci")
    parser.add_argument("metadata_out", help="opa metadata after assigning loci")
    parser.add_argument("locus_assignment_ranges", help="output file path for locus assignment ranges")
    parser.add_argument("png_filename", help="png filename for figure of opa assignment definitions")
    parser.add_argument("pdf_filename", help="pdf filename for figure of opa assignment definitions")
    return parser.parse_args()

args = get_args()

# Load data
opa_metadata = pd.read_csv(args.metadata_in, index_col = 0)
opa_metadata.sort_values(by = ['strain', 'start_reordered_flipped'], inplace = True)

# Define ranges where the opa loci will be labeled
starts = [0, 0.8*10**6, 1.2*10**6, 1.4*10**6, 1.8*10**6, 2.0*10**6]
stops = [0.15*10**6, 1.1*10**6, 1.3*10**6, 1.6*10**6, 1.9*10**6, 2.1*10**6]
num = [2, 3, 1, 3, 1, 1]
labels = ['A,B', 'C,D,E', 'F', 'G,H,I', 'J', 'K']

opa_assignments = pd.DataFrame({'range_start':starts, 'range_stop':stops, 'num_opa':num, 'labels':labels})

# Assign opa loci

# Loop through strains
for strain, df in opa_metadata.groupby('strain'):
    # print(strain, df)
    
    # Keep only opa loci that have the full gene
    df.dropna(subset = ['id'], inplace = True)

    # Loop through opa loci ranges
    for i, row in opa_assignments.iterrows():

        # Find the opa genes that fall within this range
        df_bool=(df['start_reordered_flipped']>row['range_start'])&(df['start_reordered_flipped']<row['range_stop'])
        
        # Get their indices
        idx = df_bool[df_bool].index
        
        # Get the labels for the opa loci in this range
        labels = row['labels'].split(',')
        
        # Loop through labels and assign to the opa genes in order
        for j in range(len(labels)):
            if j<len(idx):
                opa_metadata.at[idx[j], 'locus'] = labels[j]
                
# Save locus assignments
opa_metadata.to_csv(args.metadata_out)

# Save genomic ranges for locus assignments
opa_assignments.to_csv(args.locus_assignment_ranges)

# Plot the defined ranges
plt.figure(figsize = (8,3))
plt.hist(opa_metadata['start_reordered_flipped'], 200, density = True, color = 'gray')
plt.yscale('log')

ylim = plt.gca().get_ylim()
plt.fill_between([0, 0.15*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)
plt.fill_between([0.8*10**6, 1.1*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)
plt.fill_between([1.2*10**6,1.3*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)
plt.fill_between([1.4*10**6,1.6*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)
plt.fill_between([1.8*10**6,1.9*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)
plt.fill_between([2.0*10**6,2.1*10**6], ylim[0], ylim[1], color = 'gray', alpha = 0.2)

plt.text(0.04*10**6, 10**-4.7, 'A,B', fontsize = 20)
plt.text(0.8*10**6, 10**-4.7, 'C,D,E', fontsize = 20)
plt.text(1.2*10**6, 10**-4.7, 'F', fontsize = 20)
plt.text(1.4*10**6, 10**-4.7, 'G,H,I', fontsize = 20)
plt.text(1.8*10**6, 10**-4.7, 'J', fontsize = 20)
plt.text(2.0*10**6, 10**-4.7, 'K', fontsize = 20)

plt.xticks([0, 0.5*10**6, 1*10**6, 1.5*10**6, 2*10**6], ['0', '0.5 Mbp', '1 Mbp', '1.5 Mbp', '2 Mbp'], fontsize = 15)
plt.yticks(fontsize = 15)
plt.ylim(ylim)
plt.xlabel('Rearranged genomic position', fontsize = 20)
plt.ylabel('Density of\n$opa$ genes', fontsize = 20)
plt.tight_layout()
plt.savefig(args.png_filename, dpi = 300)
plt.savefig(args.pdf_filename)
plt.show()