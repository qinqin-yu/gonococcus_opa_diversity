#!/usr/bin/env python
import glob
import os
import argparse

import pandas as pd
import numpy as np

### DEFINE FUNCTIONS ###

def get_args():
    parser = argparse.ArgumentParser(description='Get position of opa genes in genome, mapped to FA1090 genome')
    parser.add_argument("xmfa_dir", help="XMFA files directory from Progressivemauve output")
    parser.add_argument("opa_locations", help="Directory of opa locations in original genome")
    parser.add_argument("outpath_lcbs", help="Output path for locally collinear blocks information")
    parser.add_argument("outpath_opa_metadata", help="Output path for opa metadata (locus, position in genome, etc)")
    return parser.parse_args()

def parse_pairwise_progressivemauve_output(xmfa_filename):
    # This function takes the xmfa file output from running progressivemauve on two genomes (reference and query) 
    # and parses it to get the location and strand of the Locally Collinear Blocks (LCB) in the reference and query
    
    # Load xmfa file
    file = open(xmfa_filename, 'r')
    lines = file.readlines()
    
    # Filter to only lines that have the relevant information (discard sequences)
    filtered_lines = [x for x in lines if '>' in x]

    strains = []
    starts = []
    stops = []
    strands = []

    # Loop through lines and get the start, stop, and strand of the LCB
    for line in filtered_lines:
        line_split = line.split(' ')
        line_split2 = line_split[1].split(':')
        strain = int(line_split2[0])
        start_stop = line_split2[1].split('-')
        strains.append(strain)
        starts.append(int(start_stop[0]))
        stops.append(int(start_stop[1]))
        strand = line_split[2]
        if strand == '+':
            strands.append(1)
        elif strand == '-':
            strands.append(-1)
    df = pd.DataFrame({'strain':strains,
                 'start':starts,
                 'stop':stops,
                 'strand':strands})
    
    # Since the xmfa file saves the blocks in the genome that don't fall into an LCB, we want to remove those.
    # The file first lists the LCB alternating between strains 1 and 2, and then lists all of the blocks from strain 1 that didn't
    # fall within an LCB followed by all of the blocks from strain 2 that didn't fall within an LCB.
    
    # Check first if there are any blocks that don't fall into an LCB
    if len(np.where(df['strain'].diff()==0)[0])>0:
        idx = np.where(df['strain'].diff()==0)[0][0]-1
        df = df.iloc[:idx]

    # Save into dataframe
    
    df_1 = df[df['strain']==1].reset_index(drop = True)
    df_1.rename({'start':'ref_lcb_start', 'stop':'ref_lcb_stop', 'strand':'query_lcb_strand'}, axis = 'columns', inplace = True)
    df_2 = df[df['strain']==2].reset_index(drop = True)
    df_2.rename({'start':'query_lcb_start', 'stop':'query_lcb_stop', 'strand':'ref_lcb_strand'}, axis = 'columns', inplace = True)

    # Note that I'm switching the reference and query lcb strands because progressivemauve seems to assign the query as always having + strand.
    
    lcbs = pd.concat([df_1, df_2], axis = 1)
    lcbs.drop(['strain'], axis = 'columns', inplace = True)
    lcbs = lcbs[['ref_lcb_start', 'ref_lcb_stop', 'ref_lcb_strand', 'query_lcb_start', 'query_lcb_stop', 'query_lcb_strand']]
    lcbs['lcb_num'] = lcbs.index

    # Calculate the ratio of the lengths of the lcb in the reference and query (used for scaling the locations of opa genes later)
    lcbs['length_ratio'] = (lcbs['ref_lcb_stop']-lcbs['ref_lcb_start'])/(lcbs['query_lcb_stop']-lcbs['query_lcb_start'])

    return lcbs

def reorder_opa_locations(opa_locations, lcbs):
    # This function takes the opa locations and reorders and flips them according to the LCBs
    
    for i, row in opa_locations.iterrows():
        
        # Match the opa locations to the LCB
        lcb_matches = np.where((lcbs['query_lcb_start']<row['start'])&(lcbs['query_lcb_stop']>row['start']))[0]
        
        # Check that there is only one match
        if len(lcb_matches)==1:
            lcb_match = int(lcb_matches[0])
            opa_locations.at[i, 'lcb_num'] = lcb_match
        elif len(lcb_matches)==0:
            print(row['id'] + ' at position ' + str(round(row['start'])) + ' is not in any LCB')
        elif len(lcb_matches)>1:
            print(row['id'] + ' at position ' + str(round(row['start'])) + ' is found in more than one LCB')

    # Merge the info about the LCB
    opa_locations = opa_locations.merge(lcbs, on = 'lcb_num', how = 'left')
    
    # Shift the opa position based on LCB
    # Do this by calculating the opa start position by shifting it to where it is relative to the LCB in the reference strain
    # (shift and rescale)
    opa_locations['start_reordered'] = round((opa_locations['start'] - opa_locations['query_lcb_start'])*opa_locations['length_ratio']+opa_locations['ref_lcb_start'])

    # Check if any of the LCBs need to be flipped (due to inversions)
    opa_locations['query_lcb_strand_multiplier'] = 0
    opa_locations.loc[opa_locations['query_lcb_strand']==-1, 'query_lcb_strand_multiplier'] = 1
    
    opa_locations['start_reordered_flipped'] = opa_locations['start_reordered'] + opa_locations['query_lcb_strand_multiplier']*(opa_locations['ref_lcb_start'] + opa_locations['ref_lcb_stop'] - 2*opa_locations['start_reordered'])
    opa_locations['strand_flipped'] = opa_locations['strand']*opa_locations['query_lcb_strand']
    
    # Fill nan values
    opa_locations.loc[np.isnan(opa_locations['start_reordered']), 'start_reordered'] = opa_locations['start']
    opa_locations.loc[np.isnan(opa_locations['start_reordered_flipped']), 'start_reordered_flipped'] = opa_locations['start']
    opa_locations.loc[np.isnan(opa_locations['strand_flipped']), 'strand_flipped'] = opa_locations['strand']
    
    opa_locations.drop(labels = ['query_strain', 'ref_strain', 'query_lcb_strand_multiplier'], axis = 'columns', inplace = True)
    
    return opa_locations

### RUN CODE ###

args = get_args()

xmfa_filenames = glob.glob(args.xmfa_dir + '*.xmfa')
opa_locations_all = pd.DataFrame()
lcbs_all = pd.DataFrame()
for xmfa_filename in xmfa_filenames:
    basename = os.path.basename(xmfa_filename)
    strain = basename[:basename.find('.xmfa')]
    if strain!='FA1090': # Do for all except reference strain
        lcbs = parse_pairwise_progressivemauve_output(xmfa_filename)
        lcbs['ref_strain'] = 'FA1090'
        lcbs['query_strain'] = strain
        
        opa_locations = pd.read_csv(args.opa_locations + strain + '.csv', index_col = 0)
        opa_locations.dropna(subset = 'id', inplace = True, ignore_index = True)
        opa_locations = reorder_opa_locations(opa_locations, lcbs)
        opa_locations['strain'] = strain
        
        lcbs_all = pd.concat([lcbs_all, lcbs])
        opa_locations_all = pd.concat([opa_locations_all, opa_locations])

# Add opa locations for reference strain (FA1090)
opa_locations = pd.read_csv(args.opa_locations + '/FA1090.csv', index_col = 0)
opa_locations['start_reordered_flipped'] = opa_locations['start']
opa_locations['strand_flipped'] = opa_locations['strand']
opa_locations['strain'] = 'FA1090'
opa_locations_all = pd.concat([opa_locations_all, opa_locations])

# Sort strains
lcbs_all.sort_values(by = ['query_strain', 'lcb_num'], inplace = True, ignore_index = True)
opa_locations_all.sort_values(by = ['strain', 'start'], inplace = True, ignore_index = True)

# Save to file
lcbs_all.to_csv(args.outpath_lcbs, index = False)
opa_locations_all.to_csv(args.outpath_opa_metadata)