#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
font = {'family' : 'Arial',
        'size'   : 12}
matplotlib.rc('font', **font)

def get_args():
    parser = argparse.ArgumentParser(description='Plot distribution of in-frame opa genes by strain')
    parser.add_argument("metadata", help="opa metadata file")
    parser.add_argument("figure_png", help="png figure filename")
    parser.add_argument("figure_pdf", help="pdf figure filename")
    return parser.parse_args()

args = get_args()

# Frame on by strain

metadata = pd.read_csv(args.metadata, index_col = 0)

strains = []
frac_in_frame = []
num_in_frame = []
for strain, df in metadata.groupby('strain'):
    strains.append(strain)
    frac_in_frame.append((df['in_frame']==1).sum()/len(df))
    num_in_frame.append((df['in_frame']==1).sum())
df_in_frame = pd.DataFrame({'strain':strains, 'num_in_frame':num_in_frame})

# Make plots

plt.figure(figsize = (4,3))
sns.countplot(df_in_frame, x="num_in_frame", color = "#ad6a6c")
plt.xlabel('Number of in-frame $opa$ genes')
plt.ylabel('Number of strains')
plt.tight_layout()
plt.savefig(args.figure_png, dpi = 300)
plt.savefig(args.figure_pdf)
