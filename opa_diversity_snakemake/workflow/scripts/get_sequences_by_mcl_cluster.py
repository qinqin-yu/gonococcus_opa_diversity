import pandas as pd
import os
from Bio import SeqIO
import mcl_functions as mf


out_path = '../data/mcl/cluster_sequences/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
    
regions = ['hv1', 'hv2', 'sv', 'con']
inflations = [1.4, 2, 4, 6, 8, 10, 12, 14]

out_filenames = []
for region in regions:
    record_dict = SeqIO.to_dict(SeqIO.parse("../data/mash/opa_sequences_no_repeats_" + region + "_ungapped.fa", "fasta"))
    for inflation in inflations:
        clusters = mf.parse_mcl_clusters('../data/mcl/dump.' + region + '.nci.I' + str(int(inflation*10)))
        for cluster, df in clusters.groupby('cluster'):
            cluster_sequences = []
            for i, row in df.iterrows():
                cluster_sequences.append(record_dict[row['id']])
            out_filename = region + '_I' + str(int(inflation*10)) + '_cluster_' + str(cluster) + '.fa'
            SeqIO.write(cluster_sequences, out_path + out_filename, "fasta")
            out_filenames.append(out_filename)

with open(out_path + 'filenames.txt', 'w') as f:
    for item in out_filenames:
        f.write(item + '\n')
f.close()