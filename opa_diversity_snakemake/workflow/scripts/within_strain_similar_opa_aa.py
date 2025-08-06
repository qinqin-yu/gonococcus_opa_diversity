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

### PLOT RESULTS 

plt.figure(figsize = (22,4))
sns.scatterplot(data=similar_opas_counts, x='strain', y='num_similar', size='num_occurrences', hue = 'num_occurrences', hue_norm = (0, 4))
plt.xticks(rotation = 90)
plt.xlabel('Isolate')
plt.ylabel('Number of similar $opa$ genes')
plt.title(str(len(strains_with_similar_opa)) + '/' + str(len(all_strains)) + ' isolates have similar $opa$ genes within the isolate')
plt.legend(title = '# occurrences', loc = 'upper right')
plt.xlim([-1, len(all_strains)])
plt.tight_layout()
plt.savefig(args.png_similar_by_strain, dpi = 300)
plt.savefig(args.pdf_similar_by_strain)
plt.show()

# Make plot showing the fraction of strains with similar opa genes as a function of the cutoff distance

similarity = np.linspace(0, 1, 200)
num_strains_similar_opa = []
for i in similarity:
    thresh = 1-i #95% similarly
    num_strains_similar_opa.append(len(np.unique(within_strain_distance[within_strain_distance['distance']<=thresh]['strain_A'])))

frac_strains_similar_opa = np.array(num_strains_similar_opa)/len(np.unique(within_strain_distance['strain_A']))

print('percent of strains with pair of opa genes with at least 95% similarity: ', 100*frac_strains_similar_opa[np.where(similarity>0.95)[0][0]])
print('number of strains with pair of opa genes with at least 95% similarity: ', np.array(num_strains_similar_opa)[np.where(similarity>0.95)[0][0]])
print('total number of strains: ', len(np.unique(within_strain_distance['strain_A'])))

# Plot (with broken axis)      
fig, (ax, ax2) = plt.subplots(1, 2, sharey = True, width_ratios = [1, 5], figsize = (4,3))

ax.plot(similarity*100, frac_strains_similar_opa*100, color = 'gray')
ax2.plot(similarity*100, frac_strains_similar_opa*100, color = 'gray')

ax.set_xlim(-3, 10)
ax2.set_xlim(75, 101)

ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
ax2.tick_params(left = False)

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d, 1+d+0.12), (-d, +d), **kwargs)
ax.plot((1-d, 1+d+0.12), (1-d, 1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the right axes
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)

ax2.set_xlabel('Similarity threshold\n(% amino acid identity)')
ax.set_ylabel('% isolates with at least 2\nsimilar $opa$ genes')
ax.set_ylim([-5, 105])
plt.tight_layout()
plt.savefig(args.png_vary_thresh, dpi = 300)
plt.savefig(args.pdf_vary_thresh)
plt.show()
    
    
