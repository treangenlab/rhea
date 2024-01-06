## function to extract sequences from a MetaFlye gfa graph and into a fasta

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Extract list of nodes from gfa')
parser.add_argument('gfa', help='Path to the metaflye gfa')
parser.add_argument('nodes', help='Path to the tsv; first column is nodes to keep')
parser.add_argument('-o', '--output', default='./extracted.fa', help='Name of the output fasta (default: extracted)')
args = parser.parse_args()

df = pd.read_csv(args.nodes, sep='\t')
node_list = list(df['node'])

with open(args.output, 'w') as fasta_file:
    with open(args.gfa, "r", encoding='UTF-8') as file:
        for line in file:
            if line.startswith("S"):  # add a node
                split_line = line[:-1].split("\t")
                node_id = split_line[1]
                sequence = split_line[2]
                if node_id in node_list:
                    fasta_file.write(">{}\n{}\n".format(node_id, sequence))
