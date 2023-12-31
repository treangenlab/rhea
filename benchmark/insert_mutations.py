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

# Rest of the code remains the same as the previous example
genome_length = len(genome_sequence)
start_pos = random.randint(0, genome_length - args.max_length)
end_pos = start_pos + random.randint(args.min_length, args.max_length)

selected_region = genome_sequence[start_pos:end_pos]
deleted_length = end_pos - start_pos
inserted_length = len(selected_region)
replacement_segment = ''.join(random.choice('ACGT') for _ in range(inserted_length))

mutated_genome = genome_sequence[:start_pos] + replacement_segment + genome_sequence[end_pos:]

# Save the mutated genome sequence to a new FASTA file with the extracted header
with open(args.output, 'w') as output_file:
    output_file.write(f'>{header}\n' + mutated_genome)

# Save alteration details in the TSV file
with open(args.output_tsv, 'w') as tsv_file:
    tsv_file.write("Start_Position\tDeleted_Length\tInserted_Length\tDeleted_Seq\tInserted_Seq\n")
    tsv_file.write(f"{start_pos}\t{deleted_length}\t{inserted_length}\t{genome_sequence[start_pos:end_pos]}\t{replacement_segment}\n")

# Save inserted and deleted sequences as separate FASTA files
with open(args.deleted_output, 'w') as deleted_file:
    deleted_file.write(f'>{header}_deleted\n{genome_sequence[start_pos:end_pos]}')
with open(args.inserted_output, 'w') as inserted_file:
    inserted_file.write(f'>{header}_inserted\n{replacement_segment}')
