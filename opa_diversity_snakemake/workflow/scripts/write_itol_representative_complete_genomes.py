#!/usr/bin/env python

import pandas as pd
import glob
import seaborn as sns
import numpy as np
import argparse

import itol_functions as fnc

def get_args():
    parser = argparse.ArgumentParser(description='Write itol for number of opa genes')
    parser.add_argument("public_assemblies_folder", help="Folder that contains the publicly availble complete genome assemblies, with file endings .fa")
    parser.add_argument("lab_assemblies_folder", help="Folder that contains the lab generated complete genome assemblies, with file endings .fa")
    parser.add_argument("itol_filename", help="itol filename")
    return parser.parse_args()

args = get_args()

public_filenames = glob.glob(args.public_assemblies_folder + "/*.fa")
public_complete = []
for filename in public_filenames:
    public_complete.append(filename.split('/')[-1].split('.fa')[0].replace("#","_"))

lab_filenames = glob.glob(args.lab_assemblies_folder + "/*.fa")
lab_complete = []
for filename in lab_filenames:
    lab_complete.append(filename.split('/')[-1].split('.fa')[0].replace("#","_"))

legend = pd.DataFrame({'wgs_id':public_complete, 'Complete genome':'Publicly available'})
legend = pd.concat([legend, pd.DataFrame({'wgs_id':lab_complete, 'Complete genome':'This study'})], ignore_index = True)

# Write itol for anatomic site

sample_name = 'wgs_id'
annotation = 'Complete genome'
output_path = args.itol_filename

unique_annotations = np.unique(legend[annotation])
color_palette = 'Set2'
colors = sns.color_palette(color_palette, len(unique_annotations)).as_hex()
colors_dict = dict(zip(unique_annotations, colors))

fnc.itol_colorstrip(legend, annotation, sample_name, colors_dict, output_path)