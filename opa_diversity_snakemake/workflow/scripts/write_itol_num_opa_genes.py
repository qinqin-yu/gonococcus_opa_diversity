#!/usr/bin/env python
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import itol_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Write itol for number of opa genes')
    parser.add_argument("opa_num", help="Filename for opa numbers")
    parser.add_argument("itol_filename", help="itol filename")
    return parser.parse_args()

args = get_args()
df = pd.read_csv(args.opa_num, index_col = 0)

df['strain'] = df['strain'].str.replace('#', '_')
df.rename({'num_opa_genes':'# opa genes'}, axis = 'columns', inplace = True)
                   
legend = df[['strain', '# opa genes']]
annotation = '# opa genes'
sample_name = 'strain'
output_path = args.itol_filename

unique_annotations = np.unique(df[annotation])
color_palette = 'mako'
colors = ["#422424", "#AF6E6E", "#DBBDBD", 'white'] #'#0D0707', '#6A3939',
# colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)