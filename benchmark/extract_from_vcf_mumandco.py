## function to convert sequences from the SURVIVOR VCF output into a fasta

import os
import argparse
import subprocess
import pandas as pd


parser = argparse.ArgumentParser(description='Extract SV information frmo SURVIVOR output VCF')
parser.add_argument('VCF_files', nargs="+", help='Path to the SURVIVOR VCF file')
parser.add_argument('-o', '--output', default='./extracted.fasta', help='Name of fasta file (default: extracted)')
args = parser.parse_args()

headers = []
sequences = []

for vcf_file in args.VCF_files:
    df = pd.read_csv(vcf_file, sep='\t')
    genus = vcf_file.split("/")[-1].split(".")[0]
    new_headers = [f"{genus}-{item}-{index + 1}" for index, item in enumerate(df['SV_type'])]

    headers = headers + new_headers
    sequences = sequences + list(df['fragments'])


with open(args.output, 'w') as fasta_file:
    # Iterate through each row
    for i in range(len(headers)):
        fasta_file.write(">{}\n{}\n".format(headers[i], sequences[i]))






