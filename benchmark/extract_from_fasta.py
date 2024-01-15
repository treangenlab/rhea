## function to extract specfied sequences from a fasta

import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Extract sequences from fasta in first col of nodes_tsv')
parser.add_argument('fasta', help='Path to full fasta')
parser.add_argument('nodes_tsv', help='Path to .tsv where the first column has desired node ids')
parser.add_argument('-o', '--output', default='./extracted.fa', help='Name of the output fasta (default: extracted)')
args = parser.parse_args()

df = pd.read_csv(args.nodes_tsv, sep='\t')
keep_list = list(df.iloc[:, 0])

records = []
set_ids = set()
for record in SeqIO.parse(args.fasta, "fasta"):
    if record.id in keep_list and record.id not in set_ids:
        records.append(record)
        set_ids.add(record.id)

SeqIO.write(records, args.output, "fasta")







