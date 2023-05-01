
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_graph(graph_path):
    sequences = []
    with open(graph_path, "r") as f:

        for line in f:
            if line.startswith("S"):
                split_line = line[:-1].split("\t")
                new_node = split_line[1].split("_")[-1]
                sequences.append(SeqRecord(Seq(split_line[2]), id="edge_{}".format(new_node)))

    SeqIO.write(sequences, "graph-nodes.fa", "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_graph", type=str,
        help="path to .gfa assembly graph")

    args = parser.parse_args()

    parse_graph(args.input_graph)
