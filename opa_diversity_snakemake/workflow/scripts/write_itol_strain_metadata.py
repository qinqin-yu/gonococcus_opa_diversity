#!/usr/bin/env python
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import itol_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Write itol for number of opa genes')
    parser.add_argument("opa_metadata", help="Opa metadata filename")
    parser.add_argument("genome_metadata", help="Genome metadata filename")
    parser.add_argument("itol_anatomic_site_filename", help="itol anatomic site filename")
    parser.add_argument("itol_gender_filename", help="itol gender filename")
    parser.add_argument("itol_sexual_behavior_filename", help="itol sexual behavior filename")
    return parser.parse_args()

# Load metadata
args = get_args()
opa = pd.read_csv(args.opa_metadata, index_col = 0)
metadata = pd.read_csv(args.genome_metadata, index_col = 0)

strains = np.unique(opa['strain'])
strains_df = pd.DataFrame({'wgs_id':strains})
merged = strains_df.merge(metadata, on = 'wgs_id', how = 'left')

merged['wgs_id'] = merged['wgs_id'].str.replace('#', '_')

merged.fillna(value='NA', inplace=True)

# Write itol for anatomic site

merged.rename({'anatomic_site':'Anatomic site'}, axis = 'columns', inplace = True)

sample_name = 'wgs_id'
annotation = 'Anatomic site'
legend = merged[[sample_name, annotation]]
output_path = args.itol_anatomic_site_filename

unique_annotations = np.unique(merged[annotation])
unique_annotations = unique_annotations[unique_annotations!='NA']
color_palette = 'Set2'
colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))
colors_dict['NA'] = '#808080'

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)

# Write itol for gender

merged.rename({'gender':'Gender'}, axis = 'columns', inplace = True)

sample_name = 'wgs_id'
annotation = 'Gender'
legend = merged[[sample_name, annotation]]
output_path = args.itol_gender_filename

unique_annotations = np.unique(merged[annotation])
unique_annotations = unique_annotations[unique_annotations!='NA']
color_palette = 'icefire'
colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))
colors_dict['NA'] = '#808080'

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)

# Write itol for sexual behavior

merged.rename({'sexual_behavior':'Sexual behavior'}, axis = 'columns', inplace = True)

sample_name = 'wgs_id'
annotation = 'Sexual behavior'
legend = merged[[sample_name, annotation]]
output_path = args.itol_sexual_behavior_filename

unique_annotations = np.unique(merged[annotation])
unique_annotations = unique_annotations[unique_annotations!='NA']
color_palette = 'Paired'
colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))
colors_dict['NA'] = '#808080'

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)