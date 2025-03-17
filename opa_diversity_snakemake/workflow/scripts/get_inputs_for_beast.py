#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = argparse.ArgumentParser(description='Get dates for BEAST analysis')
    parser.add_argument("metadata", help="Genomes metadata file")
    parser.add_argument("alignment", help="Filtered polymorphic sites from gubbins")
    parser.add_argument("alignment_with_dates", help="Filtered polymorphic sites from gubbins, only samples with dates")
    parser.add_argument("dates", help="Dates file for BEAST input")
    return parser.parse_args()

args = get_args()


# Load metadata and format dates for beast

metadata = pd.read_csv(args.metadata)
metadata_dates = metadata[['wgs_id', 'date']].copy()
metadata_dates.rename({'wgs_id':'name'}, axis = 'columns', inplace = True)

date_split = metadata_dates['date'].str.split('-', expand = True)
date_split['month'] = date_split.loc[~date_split[0].isnull(), 1].fillna('1').astype('int').astype('str')
date_split['day'] = date_split.loc[~date_split[0].isnull(), 2].fillna('01')
metadata_dates['reformatted_date'] = date_split[0] + '-' + date_split['month'] + '-' + date_split['day']

# Load alignment for beast
alignment = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))
names = list(alignment.keys())
names_df = pd.DataFrame({'name':names})

# Filter metadata to strains complete genomes and with dates
metadata_beast = names_df.merge(metadata_dates, on = 'name', how = 'left')
metadata_beast.dropna(subset = ['reformatted_date'], ignore_index = True, inplace = True)
strains_with_dates = metadata_beast['name'].values

# Filter alignment to strains with dates
records = []
for name in alignment:
    # print(opa_genes[opa].seq)
    if name in strains_with_dates:
        records.append(alignment[name])

# Save alignment
SeqIO.write(records, args.alignment_with_dates, "fasta")

# Save dates
metadata_beast[['name', 'reformatted_date']].to_csv(args.dates, sep = '\t', header = False, index = False)