#!/usr/bin/env python
import argparse
import regex
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Get list of opa genes to exclude because of frameshift mutation in orf and make new sequence file')
    parser.add_argument("opa_sequences_aa_filename", help="Concatenated fasta file of all opa aa sequences")
    parser.add_argument("filtered_aa_filename", help="Fasta file of all opa aa sequences that don't have frameshift mutations")
    parser.add_argument("excluded_opas", help="csv listing the opa sequences that are excluded because of frameshift mutations")
    return parser.parse_args()

args = get_args()

# Get list of opa genes to exclude because of frameshift mutation in orf and make new sequence file
filename = args.opa_sequences_aa_filename
opa_genes = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
opa_genes_to_exclude = []
records = []
for opa_gene in opa_genes:
    # Exclude if the final 21 codons (including stop codon) have more than 4 errors OR there is a premature stop codon.
    # Checked that all of the opa genes where this sequence was found had fewer than 4 errors.
    seq = opa_genes[opa_gene].seq
    match = regex.search(r"(?e)(RLENTRFKTHEASLGVRYRF*){e<=4}", str(seq[-21:]))
    if (not match)|(str(seq).find('*')<(len(seq)-1)):
        opa_genes_to_exclude.append(opa_gene)
    else:
        records.append(opa_genes[opa_gene])
SeqIO.write(records, args.filtered_aa_filename, "fasta")
    
output_file = args.excluded_opas
with open(output_file, "w") as outfile:
    outfile.write("opa_gene\n")
    for o in opa_genes_to_exclude:
        outfile.write(o + "\n")