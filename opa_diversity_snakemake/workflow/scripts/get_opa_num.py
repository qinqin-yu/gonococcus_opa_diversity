#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Summarize the number of opa genes per strain')
    parser.add_argument("opa_metadata", help="Filename for opa metadata")
    parser.add_argument("opa_num", help="Filename for opa numbers")
    return parser.parse_args()

args = get_args()

opa_metadata = pd.read_csv(args.opa_metadata, index_col = 0)

strains = []
num_opas = []
for value, df in opa_metadata.groupby('strain'):
    strains.append(value)
    
    # Drop rows if id is missing
    df.dropna(subset = ['id'], axis=0, inplace = True)
    
    num_opas.append(len(df))
opa_summary = pd.DataFrame({'strain':strains, 'num_opa_genes':num_opas})

opa_summary.sort_values('num_opa_genes', inplace = True)
opa_summary.reset_index(inplace = True, drop = True)
opa_summary.to_csv(args.opa_num)