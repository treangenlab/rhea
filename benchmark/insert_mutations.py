## function to randomly insert long mutations into a genome fasta

import argparse
import random
from Bio import SeqIO

# Argument parser setup
parser = argparse.ArgumentParser(description='Simulate mutation in a genome segment')
parser.add_argument('genome_file', help='Path to the genome file in FASTA format')
parser.add_argument('-o', '--output', default='mutated_genome.fasta', help='Name of the output file (default: mutated_genome.fasta)')
parser.add_argument('-min', '--min_length', type=int, default=500, help='Minimum length of mutation interval (default: 100)')
parser.add_argument('-max', '--max_length', type=int, default=2000, help='Maximum length of mutation interval (default: 1000)')
parser.add_argument('-tsv', '--output_tsv', default='mutations.tsv', help='Name of the output TSV file to track alterations (default: mutations.tsv)')
parser.add_argument('-del', '--deleted_output', default='deleted_sequences.fasta', help='Name of the output file for deleted sequences (default: deleted_sequences.fasta)')
parser.add_argument('-ins', '--inserted_output', default='inserted_sequences.fasta', help='Name of the output file for inserted sequences (default: inserted_sequences.fasta)')
args = parser.parse_args()

# Read the genome sequence from the file and extract the header
header = ""
with open(args.genome_file, 'r') as genome_file:
    for record in SeqIO.parse(genome_file, 'fasta'):
        header = record.description
        genome_sequence = str(record.seq)
        break  # Only consider the first sequence if there are multiple


# Generate subsequent random numbers that satisfy the condition
N_mutations = 10
deleted_positions = []
genome_length = len(genome_sequence)
while len(deleted_positions) < N_mutations:
    start_pos = random.randint(0, genome_length - args.max_length)
    if all(abs(start_pos - x[0]) >= args.max_length for x in deleted_positions):
        end_pos = start_pos + random.randint(args.min_length, args.max_length)
        deleted_positions.append((start_pos, end_pos))
sorted_deleted_positions = sorted(deleted_positions, key=lambda x: x[0])
sorted_deleted_positions.append((genome_length, genome_length))

# generate random sequences to insert
replacement_segments = []
for i in range(N_mutations):
    inserted_length = random.randint(args.min_length, args.max_length)
    replacement_segments.append(''.join(random.choice('ACGT') for _ in range(inserted_length)))

# track deleted
deleted_segments = []
for (start_pos, end_pos) in sorted_deleted_positions:
    deleted_segments.append(genome_sequence[start_pos:end_pos])

# generate new mutated genome
mutated_genome = genome_sequence[:sorted_deleted_positions[0][0]] + \
                 replacement_segments[0] + \
                 genome_sequence[sorted_deleted_positions[0][1]:sorted_deleted_positions[1][0]]
for i in range(1, N_mutations):
    mutated_genome = mutated_genome + replacement_segments[i] + \
                     genome_sequence[sorted_deleted_positions[i][1]:sorted_deleted_positions[i+1][0]]

# Save the mutated genome sequence to a new FASTA file with the extracted header
with open(args.output, 'w') as output_file:
    output_file.write(f'>{header}-mutated\n' + mutated_genome)

# Save inserted and deleted sequences as separate FASTA files
with open(args.deleted_output, 'w') as deleted_file:
    for i in range(N_mutations):
        if i != 0:
            deleted_file.write(f'\n')
        deleted_file.write(f'>{header}-{i}_deleted\n{deleted_segments[i]}')
with open(args.inserted_output, 'w') as inserted_file:
    for i in range(N_mutations):
        if i != 0:
            inserted_file.write(f'\n')
        inserted_file.write(f'>{header}-{i}_inserted\n{replacement_segments[i]}')

# Save alteration details in the TSV file
with open(args.output_tsv, 'w') as tsv_file:
    tsv_file.write("ID\tStart_Position\tDeleted_Length\tInserted_Length\tDeleted_Seq\tInserted_Seq\n")
    for i in range(N_mutations):
        inserted = replacement_segments[i]
        deleted = deleted_segments[i]
        tsv_file.write(f"{header}-{i}\t{sorted_deleted_positions[i][0]}\t{len(deleted)}\t{len(inserted)}\t{deleted}\t{inserted}\n")
