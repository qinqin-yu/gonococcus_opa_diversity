#!/usr/bin/env python

import pandas as pd
import numpy as np

import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Plot result of MCL clustering')
    parser.add_argument("info_filename", help="Filename of info file")
    parser.add_argument("meet_dist_filename", help="Filename of meet distance")
    parser.add_argument("summary_filename", help="Filename of output summary_file")
    return parser.parse_args()

args = get_args()
# Load summary statistics of clusters
f = open(args.info_filename, "r")
text = f.read()
lines = text.split('\n===\n')

efficiency = []
massfrac = []
areafrac = []
clusters = []
inflation = []
for line in lines:
    split_line = line.split(' ')
    inflation.append(int(split_line[4].split('/')[-1].split('.')[-1][1:])/10)
    efficiency.append(float(split_line[0][11:]))
    massfrac.append(float(split_line[1][9:]))
    areafrac.append(float(split_line[2][9:]))
    clusters.append(int(split_line[5][9:]))
    
cluster_info = pd.DataFrame({'inflation':inflation, 'num_clusters':clusters, 'efficiency':efficiency, 'mass_fraction':massfrac, 'area_fraction':areafrac})
cluster_info.sort_values('inflation', inplace = True, ignore_index = True)

# Load distances between clusters
f = open(args.meet_dist_filename, "r")
text = f.read()
lines = text.split('\n')
inflation_A = []
inflation_B = []
distance = []
num_nodes = []
for line in lines:
    split_line = line.split('\t')
    
    if len(split_line)>1:
        name_A = split_line[7].split('/')[-1].split('.')[-1]
        if name_A == 'meet':
            inflation_A.append(name_A)
        else:
            inflation_A.append(int(name_A[1:])/10)

        name_B = split_line[8].split('/')[-1].split('.')[-1]
        if name_B == 'meet':
            inflation_B.append(name_B)
        else:
            inflation_B.append(int(name_B[1:])/10)

        distance.append(int(split_line[0][2:]))
        num_nodes.append(int(split_line[3][3:]))

distance = pd.DataFrame({'inflation_A':inflation_A,
                        'inflation_B':inflation_B,
                        'distance':distance,
                        'num_nodes':num_nodes})
distance['percent_difference'] = 100*distance['distance']/(2*distance['num_nodes'])

# Merge cluster info and distances
# Get the distance between clustering done with two subsequent inflation parameters
for i in range(1, len(cluster_info)):
    inflation_A = cluster_info.iloc[i-1]['inflation']
    inflation_B = cluster_info.iloc[i]['inflation']
    row = distance[(distance['inflation_A']==inflation_A)&(distance['inflation_B']==inflation_B)]
    cluster_info.loc[i,'distance'] = row['distance'].values[0]
    cluster_info.loc[i,'num_nodes'] = row['num_nodes'].values[0]
    cluster_info.loc[i,'percent_difference'] = row['percent_difference'].values[0]
    
cluster_info.to_csv(args.summary_filename)