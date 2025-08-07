#!/usr/bin/env python
import argparse

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
font = {'family' : 'Arial',
        'size': 12}

matplotlib.rc('font', **font)

import networkx as nx

def get_args():
    parser = argparse.ArgumentParser(description='Determine the opa genes within strains that share a threshold sequence identity and plot results')
    parser.add_argument("within_strain_distance_filename", help="Filename of within-strain pairwise distances in csv (input)")
    parser.add_argument("similar_opas_filename", help="Filename of csv with information on similar opas within strains (output)")
    parser.add_argument("png_similar_by_strain", help="Filename of png of similar opas by strain (output)")
    parser.add_argument("pdf_similar_by_strain", help="Filename of pdf of similar opas by strain (output)")
    parser.add_argument("png_vary_thresh", help="Filename of png of percentage of strains with similar opas by similarity threshold (output)")
    parser.add_argument("pdf_vary_thresh", help="Filename of pdf of percentage of strains with similar opas by similarity threshold (output)")
    return parser.parse_args()

args = get_args()

# Load data
within_strain_distance = pd.read_csv(args.within_strain_distance_filename, index_col = 0)
within_strain_distance['similar'] = within_strain_distance['distance']<=0.05 #95% similarly

# Make graph of nodes (opa genes) and edges (similar opa genes) for all strains
nodes = np.unique(within_strain_distance[within_strain_distance['similar']][['id_A', 'id_B']].values.flatten())

nodes_with_attributes = []
for node in nodes:
    nodes_with_attributes.append((node, {'strain':node.split('_opa_')[0]}))
    
edges = list(within_strain_distance[within_strain_distance['similar']][['id_A', 'id_B']].itertuples(index=False, name=None))

G = nx.Graph()
G.add_nodes_from(nodes_with_attributes)
G.add_edges_from(edges)

# Get subgraphs for connected components of graph and record the number of similar opas in each group
strains = []
num_similar = []
ids = []
for connected_component in nx.connected_components(G):
    num_similar.append(len(connected_component))
    ids.append(list(connected_component))
    strains.append(list(connected_component)[0].split('_opa_')[0])
similar_opas = pd.DataFrame({'strain':strains, 'num_similar':num_similar, 'ids':ids})

similar_opas.to_csv(args.similar_opas_filename)

# Get the number of occurrences of a certain number of similar opa genes for a given strain
strains = []
num_similar = []
num_occurrences = []
for values, df in similar_opas.groupby(['strain', 'num_similar']):
    strains.append(values[0])
    num_similar.append(values[1])
    num_occurrences.append(len(df))
similar_opas_counts = pd.DataFrame({'strain':strains, 'num_similar':num_similar, 'num_occurrences':num_occurrences})

# Add strains that had no similar opa genes
all_strains = np.unique(within_strain_distance['strain_A'].values)
strains_with_similar_opa = np.unique(similar_opas['strain'])
strains_without_similar_opa = np.setdiff1d(all_strains, strains_with_similar_opa)
strains_without_similar_opa_df = pd.DataFrame({'strain':strains_without_similar_opa, 'num_similar':[0]*len(strains_without_similar_opa), 'num_occurrences':[np.nan]*len(strains_without_similar_opa)})

similar_opas_counts=pd.concat([similar_opas_counts, strains_without_similar_opa_df])
similar_opas_counts.sort_values(by = ['num_similar', 'num_occurrences'], ascending = False, inplace = True)
similar_opas_counts.reset_index(inplace = True, drop = True)
    
    
