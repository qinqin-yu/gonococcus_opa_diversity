#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
from BCBio import GFF
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Convert prokka annotations from gffs to csvs')
    parser.add_argument("gff_filename", help="GFF filename")
    parser.add_argument("csv_filename", help="Output csv filename")
    return parser.parse_args()

def prokka_gff2csv(gff_filename, csv_filename):

    in_handle = open(gff_filename)
    annotations = pd.DataFrame()
    for rec in GFF.parse(in_handle):
        for i in range(len(rec.features)):
            rec_dict = {}
            rec_dict['chromosome'] = rec.id
            qualifiers = rec.features[i].qualifiers
            rec_dict.update(qualifiers)
            location = rec.features[i].location
            rec_dict['start']=float(location.start)
            rec_dict['end']=float(location.end)
            rec_dict['strand']=location.strand
            annotations = pd.concat([annotations, pd.DataFrame(rec_dict)])
    annotations = annotations[annotations['source']=='prokka']
    annotations.dropna(axis = 'columns', how = 'all', inplace = True)
    annotations.reset_index(drop = True, inplace = True)
    annotations.drop(['ID', 'Name'], axis = 'columns', inplace = True)
    in_handle.close()
    annotations.to_csv(csv_filename)
    
def main():    
    args = get_args()
    prokka_gff2csv(args.gff_filename, args.csv_filename)
    
if __name__ == "__main__":
    main()