import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
font = {'family' : 'Arial',
        'size'   : 12}
matplotlib.rc('font', **font)

# Read output of treemmer
file = open('../data/treemmer/RTL.txt', 'r')
lines = file.readlines()

RTL = []
Nseq = []
for line in lines:
    if line.find('RTL')>=0:
        RTL.append(float(line[line.find('RTL')+6:line.find('N_seq')]))
        Nseq.append(int(line[line.find('N_seq')+7:-1]))
        
# Plot relative tree length as a function of the number of leaves
plt.plot(Nseq, RTL)
plt.gca().invert_xaxis()
plt.xlabel('Number of leaves')
plt.ylabel('Relative tree length')
plt.savefig('../figures/treemmer/RTL.png', dpi = 300)
plt.show()

# Plot metadata for trimmed strains
file = open('../data/treemmer/lab_strains_hybrid_genomes.final_tree.tre_trimmed_list_X_210', "r") 
trimmed_strains = file.read().split('\n')

for variable in variables:
    metadata_trimmed = metadata[metadata['wgs_id'].isin(trimmed_strains)]
    ax = sns.countplot(data = metadata_trimmed, y=variable, color = 'tab:blue', orient = 'v')
    for p in ax.patches:
        ax.annotate(int(p.get_width()),((p.get_x() + p.get_width()), p.get_y()), xytext=(1, -18),fontsize=12,color='k',textcoords='offset points', horizontalalignment='right')
    plt.xticks(rotation=90)
    plt.xlabel('Number of isolates')
    plt.ylabel('')
    plt.title(variable)
    plt.tight_layout()
    plt.savefig('../figures/treemmer/' + variable + '_trimmed.png', dpi = 300)
    plt.show()
    
# Get the strains for nanopore sequencing (trimmed strains minus the existing hybrid genomes)
metadata = pd.read_csv('../data/treemmer/lab_strains_hybrid_genomes_metadata.csv', index_col = 0)
file = open('../data/treemmer/lab_strains_hybrid_genomes.final_tree.tre_trimmed_list_X_210', "r") 
trimmed_strains = file.read().split('\n')
metadata['for_nanopore_sequencing'] = 'N'
metadata.loc[(metadata['wgs_id'].isin(trimmed_strains))&(metadata['assembly']=='draft'), 'for_nanopore_sequencing'] = 'Y'

# Save list of strains for sequencing
metadata[metadata['for_nanopore_sequencing']=='Y']['wgs_id'].to_csv('../data/treemmer/strains_for_nanopore_sequencing.txt', index = False, header = False)