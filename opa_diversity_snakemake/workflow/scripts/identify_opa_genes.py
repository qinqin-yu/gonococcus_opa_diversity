#!/usr/bin/env python
import glob
import argparse
import regex
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = argparse.ArgumentParser(description='ID opa genes in assemblies')
    parser.add_argument("assembly_filename", help="Assembly filename")
    parser.add_argument("outdir", help="Output directory")
    return parser.parse_args()

def get_sample_name(file_name):
    # if file_name.endswith("_contigs_filtered.fa"):
    #     sample_name = file_name.replace("_contigs_filtered.fa", "")
    if file_name.endswith("_with_plasmids.fna"):
        sample_name = file_name.replace("_with_plasmids.fna", "")
    elif file_name.endswith(".fna"):
        sample_name = file_name.replace(".fna", "")
    elif file_name.endswith(".fasta"):
        sample_name = file_name.replace(".fasta", "")
    elif file_name.endswith(".fa"):
        sample_name = file_name.replace(".fa", "")
    else:
        print(f"Warning: could not parse sample name from {file_name}")
        sample_name = file_name
    sample_name = os.path.basename(sample_name)
    return(sample_name)

def get_seq_record_dict(assembly_filename):
    seq_record_dict = SeqIO.to_dict(SeqIO.parse(assembly_filename, "fasta"))
    return seq_record_dict

def find_cr_term(assembly_filename):
    """Initial search for coding repeats ("cr") (CTCTT) and unique sequence near stop codon ("term") (TGCGCTACCGCTTCTGAT)"""
    
    cr_term_locations = []
    for rec in SeqIO.parse(assembly_filename, "fasta"):
        matches_cr_forward = regex.finditer(r"((CTCTT){3,100}){s<=2}", str(rec.seq))
        matches_cr_reverse = regex.finditer(r"((AAGAG){3,100}){s<=2}", str(rec.seq))
        matches_term_forward = regex.finditer(r"(TGCGCTACCGCTTCTGAT){s<=2}", str(rec.seq))
        matches_term_reverse = regex.finditer(r"(ATCAGAAGCGGTAGCGCA){s<=2}", str(rec.seq))
        for m in matches_cr_forward:
            cr_term_locations.append([rec.id, 'cr', 1, m.start(0), m.end(0), m.fuzzy_counts[0]])
        for m in matches_cr_reverse:
            cr_term_locations.append([rec.id, 'cr', -1, m.start(0), m.end(0), m.fuzzy_counts[0]])
        for m in matches_term_forward:
            cr_term_locations.append([rec.id, 'term', 1, m.start(0), m.end(0), m.fuzzy_counts[0]])
        for m in matches_term_reverse:
            cr_term_locations.append([rec.id, 'term', -1, m.start(0), m.end(0), m.fuzzy_counts[0]])
    cr_term_locations = pd.DataFrame(cr_term_locations, columns = ['chromosome', 'type', 'strand', 'start', 'stop', 'num_mutation'])
    return cr_term_locations

def pair_cr_term(cr_term_locations):
    """Pair cr and terms to each other, fix duplicated cr sequences"""
    
    cr_term_locations.drop(['num_mutation'], axis = 'columns', inplace = True)
    cr_term_locations = cr_term_locations.sort_values('stop')
    
    # Match cr and termination sequences that are nearby (do separately for those identified on the forward and reverse strands)
    df_forward = cr_term_locations[cr_term_locations['strand']==1]
    df_back = cr_term_locations[cr_term_locations['strand']==-1]
    
    if len(df_forward)>0:
        opas_forward = pd.merge_asof(df_forward[df_forward['type']=='cr'].rename(columns = {'stop':'stop_cr'}), df_forward[df_forward['type']=='term'].rename(columns = {'stop':'stop_term'}), left_on = 'stop_cr', right_on = 'stop_term', direction = 'forward', by = ['chromosome', 'strand'], suffixes = ['_cr', '_term'])
        opas_forward.drop(['type_cr', 'type_term'], axis = 'columns', inplace = True)
        opas_forward['approx_orf_length'] = opas_forward['stop_term']-opas_forward['stop_cr']
    else:
        opas_forward = pd.DataFrame(columns = ['chromosome', 'strand', 'start_cr', 'stop_cr', 'start_term', 'stop_term', 'approx_orf_length'])
    
    if len(df_back)>0:
        opas_back = pd.merge_asof(df_back[df_back['type']=='cr'].rename(columns = {'stop':'stop_cr'}), df_back[df_back['type']=='term'].rename(columns = {'stop':'stop_term'}), left_on = 'stop_cr', right_on = 'stop_term', direction = 'backward', by = ['chromosome', 'strand'], suffixes = ['_cr', '_term'])
        opas_back.drop(['type_cr', 'type_term'], axis = 'columns', inplace = True)
        opas_back['approx_orf_length'] = opas_back['start_cr']-opas_back['start_term']
    else:
        opas_back = pd.DataFrame(columns = ['chromosome', 'strand', 'start_cr', 'stop_cr', 'start_term', 'stop_term', 'approx_orf_length'])

    # Keep the pairs of cr and termination sequences if # of nucleotides between them is <1200 bp
    opa_locs = pd.concat([opas_forward, opas_back]).reset_index(drop = True)
    opa_locs = opa_locs[opa_locs['approx_orf_length']<=1200].reset_index(drop = True)
        
    # Fix duplicates (due to additional variations in cr sequence where the cr searching code split them into two cr sequences) 
    opa_locs_fix_duplicates = pd.DataFrame(columns = ['chromosome', 'strand', 'start_cr', 'stop_cr', 'start_term', 'stop_term', 'approx_orf_length'])
    for value, df_i in opa_locs.groupby(['start_term', 'stop_term']):
        if len(df_i)>1:
            df_new = df_i.tail(1).copy()
            df_new.loc[:,'start_cr'] = df_i['start_cr'].values[0]
            df_new.loc[:,'stop_cr'] = df_i['stop_cr'].values[-1]
            opa_locs_fix_duplicates = pd.concat([opa_locs_fix_duplicates, df_new])
        else:
            opa_locs_fix_duplicates = pd.concat([opa_locs_fix_duplicates, df_i])

    # Add back the termination sequences that are alone (did not find matching repeat)
    merged_all = cr_term_locations[cr_term_locations['type']=='term'].merge(opa_locs_fix_duplicates, left_on = ['chromosome', 'start', 'strand'], right_on = ['chromosome', 'start_term', 'strand'], how = 'left', indicator = True)
    term_alone = merged_all[merged_all['_merge']=='left_only'].copy()
    
    term_alone.drop(['type', 'start_term', 'stop_term', '_merge'], axis = 'columns', inplace = True)
    term_alone.rename({'start':'start_term', 'stop':'stop_term'}, axis = 'columns', inplace = True)
    opa_locs_fix_duplicates = pd.concat([opa_locs_fix_duplicates, term_alone]).reset_index(drop = True)
    opa_locs_fix_duplicates.drop(['approx_orf_length'], axis = 'columns', inplace = True)

    opa_locs_fix_duplicates['start_cr'] = opa_locs_fix_duplicates['start_cr'].astype('float')
    opa_locs_fix_duplicates['stop_cr'] = opa_locs_fix_duplicates['stop_cr'].astype('float')
    opa_locs_fix_duplicates['start_term'] = opa_locs_fix_duplicates['start_term'].astype('float')
    opa_locs_fix_duplicates['stop_term'] = opa_locs_fix_duplicates['stop_term'].astype('float')
    
    return opa_locs_fix_duplicates

def find_missing_cr(opa_metadata, seq_record_dict):
    """For termination seaquences that are alone, do a more lenient search in the area of the genome that I expect it to be in"""
    
    opa_metadata_missing_cr = opa_metadata[np.isnan(opa_metadata['start_cr'].values)]
    
    if len(opa_metadata_missing_cr)>0:
        print(len(opa_metadata_missing_cr), ' opa genes with term sequence but missing cr')

        for i, row in opa_metadata_missing_cr.iterrows():

            seq = seq_record_dict[row['chromosome']].seq

            cr_locations = []
            if row['strand']==1:
                approximate_orf = str(seq[int(row['start_term']-1200):int(row['stop_term'])])
                matches = regex.finditer(r"(((CTCTT){e<=1}){3,100})", approximate_orf)
                for m in matches:
                    cr_locations.append([m.start()+row['start_term']-1200, m.end()+row['start_term']-1200])
                    
                # If still not repeats found, then look for the sequence A(x5-7)CCTT
                if len(cr_locations)==0:
                    match = regex.search(r"(A){5,7}(CCTT){e<=1}", approximate_orf)
                    if match is not None:
                        cr_locations.append([match.end()+row['start_term']-1200, match.end()+row['start_term']-1200])
            else:
                approximate_orf = str(seq[int(row['start_term']):int(row['stop_term']+1200)])
                matches = regex.finditer(r"(((AAGAG){e<=1}){3,100})", approximate_orf)
                for m in matches:
                    cr_locations.append([m.start()+row['start_term'], m.end()+row['start_term']])
                if len(cr_locations)==0:
                    match = regex.search(r"(AAGG){e<=1}(T){5,7}", approximate_orf)
                    if match is not None:
                        cr_locations.append([match.start()+row['start_term'], match.start()+row['start_term']])
            
            # Update metadata table if cr is found
            if len(cr_locations)>0:
                start_cr = cr_locations[0][0]
                stop_cr = cr_locations[-1][-1]
                opa_metadata.loc[i,'start_cr'] = start_cr
                opa_metadata.loc[i,'stop_cr'] = stop_cr
    
            print('After more lenient search, ', len(opa_metadata[np.isnan(opa_metadata['start_cr'])]), ' opa genes with term sequence beginning portion of gene')
    else:
        print('All term sequences have an associated cr or cr-upstream sequence that was found')
        
    return opa_metadata

def get_opa_locs(opa_metadata, seq_record_dict, strain):
    
    """Get the ORF for each opa gene and assign name, refine cr location, find N-terminus"""

    opa_metadata.sort_values(by = 'start_term', inplace = True)
    opa_metadata.reset_index(inplace = True, drop = True)
    
    opa_metadata[['start', 'stop', 'n_terminus', 'id']] = np.nan # First create columns (in case there were no opas found)

    index=1
    for i, row in opa_metadata.iterrows():
        seq = seq_record_dict[row['chromosome']].seq
            
        if ~np.isnan(row['start_cr']):
            
            # Get the ORF
            if row['strand']==1:
                start_idx = seq[int(row['start_cr'])-50:int(row['start_cr'])].find('ATG')
                if start_idx >= 0:
                    start = int(start_idx-50+row['start_cr'])
                else:
                    start = np.nan
                    print('start codon not found for ' + strain + str(row['stop_term']))
                stop = int(row['stop_term']-1)
            else:
                stop_idx = seq[int(row['stop_cr']):int(row['stop_cr'])+50].find('CAT')
                if stop_idx>=1:
                    stop = int(stop_idx+row['stop_cr'])+3
                else:
                    stop = np.nan
                    print('start codon not found for ' + strain + str(row['start_term']))
                start = int(row['start_term']+1)

            # Update metadata file
            opa_metadata.loc[i,'start']=start
            opa_metadata.loc[i,'stop']=stop
                
            # Assign an ID if both start and stop have been found
            if (not np.isnan(start))&(not np.isnan(stop)):
                opa_id = strain + '_opa_' + str(index)
                opa_metadata.loc[i,'id']=opa_id
                index+=1

            # Also update the location of the coding repeat sequence to be more precise

            # If start codon was found then use the orf, otherwise then approximate the orf using slightly upstrem of start of cr
            if row['strand']==1:
                if ~np.isnan(start):
                    opa_gene_seq = seq[start:stop]
                else:
                    opa_gene_seq = seq[int(row['start_cr']-40):stop]
                    start = row['start_cr']-40
            elif row['strand']==-1:
                if ~np.isnan(stop):
                    opa_gene_seq = seq[start:stop].reverse_complement()
                else:
                    opa_gene_seq = seq[start:int(row['stop_cr']+40)].reverse_complement()
                    stop = row['stop_cr']+40

            # Check for start of cr
            
            start_match = regex.search(r"(A){5,7}(CCTT){e<=1}", str(opa_gene_seq))
            if start_match is None:
                # If this sequence is not found, then the approximate start_cr estimated from earlier is used
                opa_metadata.loc[i, 'start_cr']=row['start_cr']
                start_cr = 0 # Defining for using later in script
            else:
                start_cr = start_match.end()
                
                if row['strand']==1:
                    opa_metadata.loc[i, 'start_cr']=start+start_cr
                else:
                    opa_metadata.loc[i, 'stop_cr']=stop-start_cr

            # Check for end of cr
            cr_match = regex.search(r"(((CTCTT){e<=1}){1,100})", str(opa_gene_seq[start_cr:]))
            if cr_match is None:
                # If there are no repeats, then the stop_cr is given by the start_cr
                opa_metadata.loc[i, 'stop_cr']=opa_metadata.loc[i, 'start_cr']
            else:
                approx_stop_cr = start_cr + cr_match.end()
                start_cds_match = regex.search(r"C(?e)(CG){e<=1}",str(opa_gene_seq[approx_stop_cr:]))
                stop_cr = approx_stop_cr + start_cds_match.start()

                if row['strand']==1:
                    opa_metadata.loc[i, 'stop_cr']=start+stop_cr
                else:
                    opa_metadata.loc[i, 'start_cr']=stop-stop_cr
            
            # Find N-terminus
            nterm_match = regex.search(r"(GCAAGTGA){s<=2}", str(opa_gene_seq))
            if nterm_match is None:
                n_terminus = np.nan
                opa_metadata.loc[i, 'n_terminus'] = np.nan
            else:
                n_terminus = nterm_match.start()
                if row['strand']==1:
                    opa_metadata.loc[i, 'n_terminus'] = start+n_terminus
                else:
                    opa_metadata.loc[i, 'n_terminus'] = stop-n_terminus
            
            # Unassign the ID if N-terminus has not been found
            if np.isnan(n_terminus):
                opa_metadata.loc[i,'id']=np.nan
                index-=1
        
    return opa_metadata

def frame_on(opa_metadata):

    """Check if the N-terminus is in frame wrt start codon"""
    
    opa_metadata['in_frame'] = np.nan # First create column, in case there were no opas found

    for strand, df in opa_metadata.groupby('strand'):
        if strand == 1:
            opa_metadata.loc[df.index, 'in_frame'] = (df['n_terminus']-df['start'])%3==0
        else:
            opa_metadata.loc[df.index, 'in_frame'] = (df['stop']-df['n_terminus'])%3==0

    opa_metadata.loc[opa_metadata['in_frame']==True, 'in_frame'] = 1
    opa_metadata.loc[opa_metadata['in_frame']==False, 'in_frame'] = 0
    opa_metadata.loc[opa_metadata['id'].isnull(), 'in_frame'] = np.nan

    return opa_metadata

### Save to fasta file, one file per strain

def opa_seqs_to_fasta(opa_metadata, seq_record_dict, sample_name, outdir):
    """Save opa sequences to fasta file"""
    records = []
    for i, row in opa_metadata.iterrows():
        if not pd.isna(row['id']):
            seq = seq_record_dict[row['chromosome']].seq
            if row['strand']==1:
                records.append(SeqRecord(seq[int(row['start']):int(row['stop'])], id = row['id']))
            else:
                records.append(SeqRecord(seq[int(row['start']):int(row['stop'])].reverse_complement(), id = row['id']))
    SeqIO.write(records, outdir + '/opa_sequences/'+sample_name+".fa", "fasta")

def opa_seqs_no_cr(opa_metadata, seq_record_dict, sample_name, outdir):
    records = []
    for i, row in opa_metadata.iterrows():
        if not pd.isna(row['id']):
            seq = seq_record_dict[row['chromosome']].seq
            if row['strand']==1:
                records.append(SeqRecord(seq[int(row['stop_cr']+2):int(row['stop'])], id = row['id']))
            else:
                records.append(SeqRecord(seq[int(row['start']):int(row['start_cr']-2)].reverse_complement(), id = row['id']))
    SeqIO.write(records, outdir + '/opa_sequences_no_repeats/'+sample_name+'.fa', 'fasta')
    
def save_metadata(opa_metadata, sample_name, outdir):
    """Save metadata file"""
    opa_metadata.to_csv(outdir + '/opa_locations/' + sample_name + '.csv')
    
def main():    
    args = get_args()
    print(args.assembly_filename)
    print(os.getcwd())
    sample_name = get_sample_name(args.assembly_filename)
    cr_term_locations = find_cr_term(args.assembly_filename)
    seq_record_dict = get_seq_record_dict(args.assembly_filename)

    opa_metadata = pair_cr_term(cr_term_locations)
    opa_metadata = find_missing_cr(opa_metadata, seq_record_dict)
    opa_metadata = get_opa_locs(opa_metadata, seq_record_dict, sample_name)
    opa_metadata = frame_on(opa_metadata)
    opa_seqs_to_fasta(opa_metadata, seq_record_dict, sample_name, args.outdir)
    opa_seqs_no_cr(opa_metadata, seq_record_dict, sample_name, args.outdir)
    save_metadata(opa_metadata, sample_name, args.outdir)
    
if __name__ == "__main__":
    main()

