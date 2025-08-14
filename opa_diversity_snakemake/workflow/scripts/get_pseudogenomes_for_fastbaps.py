#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_args():
    parser = argparse.ArgumentParser(description='Get pseudogenome alignment (and replace # with _')
    parser.add_argument("pseudogenome_path", help="List of pseudogenome paths")
    parser.add_argument("pseudogenome_alignment", help="Filename of pseudogenome alignment")
    return parser.parse_args()

args = get_args()

new_records = []
for record in SeqIO.parse(args.pseudogenome_path, "fasta"):
    new_records.append(SeqRecord(
        Seq(record.seq),
        id=record.id.replace('#','_'), # to match with the gubbins tree having replace '#' by '_'
        name=record.name,
        description=record.description))
        
SeqIO.write(new_records, args.pseudogenome_alignment, "fasta")