#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import hv1_hv2_association_functions as assoc_fnc
from munkres import Munkres, print_matrix

def get_args():
    parser = argparse.ArgumentParser(description='HV1 and HV2 association when subsampling 1 genome per baps cluster')
    parser.add_argument("fastbaps_clusters", help="Fastbaps clusters")
    parser.add_argument("mcl_clusters", help="MCL clusters")
    parser.add_argument("outfile", help='Output filename')
    return parser.parse_args()

# Determines the overlap of HV1 and HV2 types, similar to 'hv1_hv2_association', but only sampling
# one strain per baps cluster to account for population structure

args = get_args()

# Import whole genome baps clusters
baps_clusters = pd.read_csv(args.fastbaps_clusters)

# Import MCL opa variable region clusters
mcl_clusters = pd.read_csv(args.mcl_clusters, index_col = 0)

# Drop reference opa
mcl_clusters=mcl_clusters[mcl_clusters['id']!='FA1090_opa_1_reference']
mcl_clusters.reset_index(inplace = True, drop = True)

strain = []
opa_id_num = []
for i in mcl_clusters['id']:
    parsed = assoc_fnc.parse_opa_name(i)
    strain.append(parsed[0])
    opa_id_num.append(int(parsed[1]))
    
mcl_clusters['strain'] = strain
mcl_clusters['opa_id_num'] = opa_id_num

locus_1_name = 'hv1_cluster'
locus_2_name = 'hv2_cluster'

# Drop identical opas in the same strain (same SV, HV1, HV2)
# The motivation is that identical opas may indicate recent duplication events which would not
# yet have had enough time to be subject to selection
mcl_clusters = mcl_clusters.drop_duplicates(subset = ['strain', 'sv_cluster', 'hv1_cluster', 'hv2_cluster'], ignore_index = True)

# Compare the actual and randomized trace when drawing one strain per baps cluster
num_draws = 100
frac_higher_trace = []

for i in range(num_draws):
    # Choose one strain per baps cluster
    strains_random_order = assoc_fnc.get_strains_random_order(baps_clusters)
    clusters = mcl_clusters[mcl_clusters['strain'].isin(strains_random_order)]

    # Calculate the actual matrix of HV1 and HV2 associations and maximize trace
    pivoted_rearranged, total = assoc_fnc.clusters_to_trace(clusters, locus_1_name, locus_2_name)

    # Randomize HV2 clusters for each opa gene and repeat above calculation
    total_randomized_all_trials = assoc_fnc.get_randomized_matrix_traces(clusters, locus_1_name, locus_2_name)
    
    # Calculate the percentage of random matrix traces that are higher than the actual trace
    frac_higher_trace.append((total_randomized_all_trials>total).sum()/len(total_randomized_all_trials))
    
frac_higher_trace = np.array(frac_higher_trace)

# Print output

num_draws_frac_higher_trace = len(frac_higher_trace[frac_higher_trace>0])

with open(args.outfile, "w") as f:
    f.write('Out of ' + str(num_draws) + ' random draws of strains (1 per baps cluster), in ' + str(num_draws_frac_higher_trace) + ' draws, the random trace was higher than the actual value.')
    if num_draws_frac_higher_trace>0:
        f.write('In each of those draws, the fraction of trials with higher trace: ', frac_higher_trace[frac_higher_trace>0])