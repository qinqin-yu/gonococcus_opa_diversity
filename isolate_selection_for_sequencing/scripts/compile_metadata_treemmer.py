import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
font = {'family' : 'Arial',
        'size'   : 12}
matplotlib.rc('font', **font)

# Read in metadata
metadata = pd.read_csv('../data/treemmer/lab_strains_hybrid_genomes_metadata.csv', index_col = 0)
metadata['year'] = metadata['date'].str.split('-', expand = True)[0]
metadata['decade'] = metadata['year'].astype('float') - (metadata['year'].astype('float')%10)

# Variables to plot
variables = ['continent', 'country', 'decade', 'gender', 'sexual_behavior', 'anatomic_site', 'assembly']

# Make plots of the number of isolates with each metadata category in the lab strains and hybrid genomes
for variable in variables:
    metadata_annotated = metadata.copy()
    ax = sns.countplot(data = metadata, y=variable, color = 'tab:blue', orient = 'v')
    for p in ax.patches:
        ax.annotate(int(p.get_width()),((p.get_x() + p.get_width()), p.get_y()), xytext=(1, -18),fontsize=12,color='k',textcoords='offset points', horizontalalignment='right')
    plt.xticks(rotation=90)
    plt.xlabel('Number of isolates')
    plt.ylabel('')
    plt.title(variable)
    plt.tight_layout()
    plt.savefig('../figures/treemmer/' + variable + '.png', dpi = 300)
    plt.show()
    
# Save the metadata in the format required by treemmer

metadata_treemmer = metadata.copy()
metadata_treemmer = metadata_treemmer[metadata_treemmer['wgs_id']!='UMASS-DGI_13']
metadata_treemmer = metadata_treemmer[metadata_treemmer['wgs_id']!='F62']
metadata_treemmer.loc[metadata_treemmer['wgs_id']=='28Bl', 'wgs_id'] = '28B'

wgs_ids = np.array([])
metadata_labels = np.array([])
for variable in variables:
    wgs_ids = np.concatenate([wgs_ids, metadata_treemmer['wgs_id'].values])
    metadata_labels = np.concatenate([metadata_labels, metadata_treemmer[variable].values])
metadata_treemmer_input = pd.DataFrame({'wgs_id':wgs_ids, 'metadata_labels':metadata_labels})

metadata_treemmer_input.to_csv('../data/treemmer/treemmer_metadata_input.csv', header = None, index = None)