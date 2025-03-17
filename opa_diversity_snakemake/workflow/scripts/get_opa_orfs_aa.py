#!/usr/bin/env python
import os
import glob
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
    
def get_args():
    parser = argparse.ArgumentParser(description='Translate opa sequences')
    parser.add_argument("nt_filename", help="Fasta file of opa nucleotide sequences without repeats")
    parser.add_argument("aa_filename", help="Fasta file of translated opa sequences")
    return parser.parse_args()

args = get_args()

filename = args.nt_filename
strain = filename.split('/')[-1][:-3]
records = []
opa_genes = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
for opa in opa_genes:
    # print(opa_genes[opa].seq)
    aa = opa_genes[opa].seq.translate()
    records.append(SeqRecord(aa, id = opa))
SeqIO.write(records, args.aa_filename, "fasta")