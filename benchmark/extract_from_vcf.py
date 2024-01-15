## function to convert sequences from the SURVIVOR VCF output into a fasta

import os
import argparse
import subprocess
import pandas as pd


parser = argparse.ArgumentParser(description='Extract SV information from SURVIVOR output VCF')
parser.add_argument('reference', help='Path to the reference genome file in FASTA format')
parser.add_argument('VCF_file', help='Path to the SURVIVOR VCF file')
parser.add_argument('-o', '--output-dir', default='./extracted', help='Name of the output directory (default: extracted)')
args = parser.parse_args()

HEADERS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with open(args.VCF_file, 'r') as file:
    lines = [line.strip() for line in file if not line.startswith('#')]

df = pd.DataFrame([line.split('\t') for line in lines], columns=HEADERS)
df['POSEND'] = df['INFO'].str.extract(r'END=(.*?);')


out_base = os.path.splitext(os.path.basename(args.reference))[0].split("_")[0]
fasta_outpath = "{}/{}-SVs.fasta".format(args.output_dir, out_base)
df_outpath = "{}/{}-SVS.tsv".format(args.output_dir, out_base)

with open(fasta_outpath, 'w') as fasta_file:
    # Iterate through each row
    for index, row in df.iterrows():
        if row['ID'][:3] == "INS":
            extracted_seq = row['ALT']
            fasta_file.write(">{}-{}\n{}\n".format(out_base, row['ID'], extracted_seq))
        elif row['ID'][:3] == "DEL":
            extracted_seq = row['REF']
            fasta_file.write(">{}-{}\n{}\n".format(out_base, row['ID'], extracted_seq))
        elif row['ID'][:3] == "DUP" or row['ID'][:3] == "INV":
            extracted_seq = subprocess.check_output("sed '1d' {} | tr -d '\n' | cut -c {}-{}"
                                                .format(args.reference, row['POS'], row['POSEND']),
                                                shell=True).decode('utf-8').split('\n')[0]
            df.at[index, 'REF'] = extracted_seq
            fasta_file.write(">{}-{}\n{}\n".format(out_base, row['ID'], extracted_seq))
df.to_csv(df_outpath, sep='\t', index=False)






