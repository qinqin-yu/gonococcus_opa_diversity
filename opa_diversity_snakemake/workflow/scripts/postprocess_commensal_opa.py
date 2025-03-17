#!/usr/bin/env python

import pandas as pd
import glob
import os

import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Combine opa locations for commensals and merge with species information')
    parser.add_argument("commensal_metadata", help="Metadata of commensal genomes")
    parser.add_argument("opa_locations_folder", help="opa locations folder")
    parser.add_argument("opa_commensals_metadata", help="opa commensals output file with speces info")
    return parser.parse_args()

args = get_args()

# Load metadata
metadata = pd.read_csv(args.commensal_metadata, index_col = 0)
metadata = metadata[['strain', 'species']]

# Replace weird characters in strain names to match what they were for the fasta filenames
metadata['strain'] = metadata['strain'].str.replace('/', '_')
metadata['strain'] = metadata['strain'].str.replace(' ', '')

# Load in all opa locations and combine into one dataframe
filenames = glob.glob(args.opa_locations_folder + '/*.csv')
opa = pd.DataFrame()
for filename in filenames:
    strain = os.path.basename(filename)[:-4]
    df = pd.read_csv(filename, index_col = 0)
    df['strain'] = strain
    opa = pd.concat([opa, df])

# Merge with metadata to get species information
opa = opa.merge(metadata, on = 'strain')

# Save results
opa.to_csv(args.opa_commensals_metadata)