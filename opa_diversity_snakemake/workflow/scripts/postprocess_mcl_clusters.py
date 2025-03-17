#!/usr/bin/env python

import mcl_functions as mf
import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Get final clusters from MCL')
    parser.add_argument("sv_filename", help="Filename of SV clusters with chosen inflation value from MCL")
    parser.add_argument("hv1_filename", help="Filename of HV1 clusters with chosen inflation value from MCL")
    parser.add_argument("hv2_filename", help="Filename of HV2 clusters with chosen inflation value from MCL")
    parser.add_argument("clusters_filename", help="Output filename of MCL clusters")
    return parser.parse_args()

args = get_args()

clusters_filenames = [args.sv_filename,
                      args.hv1_filename, 
                      args.hv2_filename]
i = 0
for clusters_filename in clusters_filenames:
    df = mf.parse_mcl_clusters(clusters_filename)
    print(clusters_filename)
    region = clusters_filename.split('.')[-3]
    df.rename({'cluster':region + '_cluster'}, axis = 1, inplace = True)
    if i==0:
        clusters = df.copy()
    else:
        clusters = clusters.merge(df, on = 'id')
    i+=1
clusters.sort_values('id', inplace = True)
clusters.reset_index(inplace = True, drop = True)
clusters.to_csv(args.clusters_filename)