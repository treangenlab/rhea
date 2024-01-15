## Script to get all the sequence ids from a fasta into a .tsv file
## used to get sequences for specific nodes from a graph fasta

import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='extract reference sequence ids')
parser.add_argument('fasta', help='Path FASTA sequences')
parser.add_argument('-o', '--output', default='./node-ids.tsv', help='Name of the output .tsv')
args = parser.parse_args()

ids = []
for record in SeqIO.parse(args.fasta, "fasta"):
    id_split = record.id.split("_")
    new_id = f"{id_split[1]}_{id_split[2][:-1]}"
    ids.append(new_id)
df = pd.DataFrame(ids)
df.to_csv(args.output, sep='\t', index=False, header=False)