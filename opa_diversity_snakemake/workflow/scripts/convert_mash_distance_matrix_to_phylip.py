#!/usr/bin/env python

import argparse

import pandas as pd
import numpy as np
    
def get_args():
    parser = argparse.ArgumentParser(description='Convert mash distance matrix to phylip format')
    parser.add_argument("mash_distance", help="Mash distance filename")
    parser.add_argument("phylip", help="Phylip filename")
    return parser.parse_args()

args = get_args()

df = pd.read_csv(args.mash_distance, delimiter = '\t', header = None, names = ['reference_id', 'query_id', 'mash_distance', 'pvalue', 'matching_hashes'])

n = int(np.sqrt(len(df)))

ids = df['reference_id'][0:n].values

phylip_filename = args.phylip
with open(phylip_filename, 'w') as outfile:
    outfile.write('\t' + str(int(n)) + '\n')
    for query_id in ids:
        mash_distances = df[df['query_id']==query_id]['mash_distance'].values
        outfile.write(query_id + '\t' + '\t'.join(mash_distances.astype('str')) + '\n')