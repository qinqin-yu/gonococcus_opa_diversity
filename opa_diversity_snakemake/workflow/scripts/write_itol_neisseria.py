#!/usr/bin/env python
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import itol_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Write itol for Neisseria species for tree of all Neisseria opa')
    parser.add_argument("opa_commensals_metadata", help="Commensal Neisseria opa metadata filename")
    parser.add_argument("opa_gc_metadata", help="GC opa metadata filename")
    parser.add_argument("itol", help="itol filename")
    return parser.parse_args()

args = get_args()

# Load metadata

metadata_commensals = pd.read_csv(args.opa_commensals_metadata, index_col = 0)
metadata_commensals.dropna(subset = ['id'], inplace = True, ignore_index = True)
metadata_commensals = metadata_commensals[['id', 'species']]

metadata_gc = pd.read_csv(args.opa_gc_metadata, index_col = 0)
metadata_gc.dropna(subset = ['id'], inplace = True, ignore_index = True)
metadata_gc = metadata_gc[['id']]
metadata_gc['species'] = 'Neisseria gonorrhoeae'

metadata = pd.concat([metadata_commensals, metadata_gc], ignore_index = True)

# Write itol for species

sample_name = 'id'
annotation = 'species'
legend = metadata[[sample_name, annotation]]
output_path = args.itol

unique_annotations = np.unique(metadata[annotation])
color_palette = 'husl'
colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)