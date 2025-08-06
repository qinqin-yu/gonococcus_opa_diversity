#!/usr/bin/env python

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import numpy as np

import os
import argparse

# Gets the different regions of the alignment
# Defined in FA1090_opa_1 by comparing the nucleotide sequence to that from Bhat et al.
# The regions in the other opa sequences are then defined using the alignment

# sv = semivariable region
# hv1 = hypervariable region 1
# hv2 = hypervariable region 2
# con = conserved region (everything else)

def get_args():
    parser = argparse.ArgumentParser(description='Gets the different regions of the alignment')
    parser.add_argument("alignment_filename", help="Nucleotide alignment (without repeats) filename")
    parser.add_argument("sv_filename", help="Output filename of semivariable region alignment")
    parser.add_argument("hv1_filename", help="Output filename of hypervariable region 1 alignment")
    parser.add_argument("hv2_filename", help="Output filename of hypervaraible region 2 alignment")
    parser.add_argument("sv_filename_ungapped", help="Output filename of semivariable region alignment with gaps removed")
    parser.add_argument("hv1_filename_ungapped", help="Output filename of hypervariable region 1 alignment with gaps removed")
    parser.add_argument("hv2_filename_ungapped", help="Output filename of hypervaraible region 2 alignment with gaps removed")
    return parser.parse_args()

args = get_args()

alignment = AlignIO.read(args.alignment_filename, 'fasta')

# Get the reference sequence
ref = str(alignment[-1].seq)

# Get the positions of in the alignment that correspond to the different variable regions
# Replace all bases by 's'
ref_replaced = ref.replace('a', 's')
ref_replaced = ref_replaced.replace('t', 's')
ref_replaced = ref_replaced.replace('c', 's')
ref_replaced = ref_replaced.replace('g', 's')

# Count number of bases (excluding gaps)
count = 0
cumcount = []
for i in ref_replaced:
    if i == 's':
        count+=1
    cumcount.append(count)
cumcount = np.array(cumcount)

# Get start and end of each variable region in alignment (including gaps)
sv_start = np.where(cumcount>=93)[0][0]
sv_end = np.where(cumcount>=122+1)[0][0]

hv1_start = np.where(cumcount>=264)[0][0]
hv1_end = np.where(cumcount>=336+1)[0][0]

hv2_start = np.where(cumcount>=456)[0][0]
hv2_end = np.where(cumcount>=598+1)[0][0]

# Get alignments and save them

sv_alignment = alignment[:, sv_start:sv_end]
hv1_alignment = alignment[:, hv1_start:hv1_end]
hv2_alignment = alignment[:, hv2_start:hv2_end]

AlignIO.write(sv_alignment, args.sv_filename, 'fasta')
AlignIO.write(hv1_alignment, args.hv1_filename, 'fasta')
AlignIO.write(hv2_alignment, args.hv2_filename, 'fasta')

# Remove gaps

records_ungapped = []
for record in sv_alignment:
    records_ungapped.append(SeqRecord(Seq(record.seq.replace("-","")),id=record.id))
SeqIO.write(records_ungapped, args.sv_filename_ungapped, "fasta")

records_ungapped = []
for record in hv1_alignment:
    records_ungapped.append(SeqRecord(Seq(record.seq.replace("-","")),id=record.id))
SeqIO.write(records_ungapped, args.hv1_filename_ungapped, "fasta")

records_ungapped = []
for record in hv2_alignment:
    records_ungapped.append(SeqRecord(Seq(record.seq.replace("-","")),id=record.id))
SeqIO.write(records_ungapped, args.hv2_filename_ungapped, "fasta")