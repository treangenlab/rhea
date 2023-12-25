#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""

import os
import logging
import argparse
import subprocess
import pandas as pd

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'Jan 2024'

PAF_headers = ['query', 'query length', 'query start', 'query end', 'relative', 'target name',
               'target length', 'target start', 'target end', 'n matches', 'alignment length',
               'quality', 'tp', 'cm', 's1', 's2', 'dv']


def count_bps(sequence_file):
    """
    Parses .gfa file from MetaFlye output into a networkx graph. Nodes are weighted by length of
    sequence. edges are weighted by coverage.

    :str graph_path: path to .gfa MetaFlye output
    :return: Weighted Networkx graph (G), dictionary of edge id to SeqRecord sequence (sequences)
    """
    base_file_ext = os.path.splitext(sequence_file.split(".gz")[0])[1]
    regex_command = '/^>/ {{next}}' if base_file_ext[-1] == "a" else 'NR%4==2'
    print_command = 'gzcat' if sequence_file.endswith('.gz') else 'cat'
    try:
        with subprocess.Popen([print_command, sequence_file], stdout=subprocess.PIPE) \
                as zcat_process:
            with subprocess.Popen(['awk', '{} {{bp += length($0)}} END {{print bp}}'
                    .format(regex_command)],
                              stdin=zcat_process.stdout, stdout=subprocess.PIPE, text=True) \
                    as awk_process:
                zcat_process.stdout.close()  # Close the unused end of the pipe
                process_output, _ = awk_process.communicate()
                return int(process_output.strip())
    except subprocess.CalledProcessError:
        print("Unable to gather bp count from input reads. No weights will be applied.")
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    parser.add_argument(
        "input", type=str, nargs='+',
        help="path to metagenome sequences files to compare, in order")
    parser.add_argument(
        "--input-graph", type=str,
        help="path to .gfa assembly graph by MetaFlye")
    parser.add_argument(
        "--bp-table", type=str,
        help="path to .tsv with bps for each input read")
    parser.add_argument(
        '--type', '-x', choices=['pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
                                 'nano-raw', 'nano-corr', 'nano-hq'],
        default='nano-raw', help='specify type for flye [nano-raw]')
    parser.add_argument(
        "--flye-exec", type=str, default='flye',
        help="path to flye executable")
    parser.add_argument(
        "--minigraph-exec", type=str, default='minigraph',
        help="path to minigraph executable")
    parser.add_argument(
        '--output-dir', '-o', type=str, default=os.path.join(os.getcwd(), "RHEA_results"),
        help='output directory name [./RHEA_results]')
    parser.add_argument(
        '--threads', type=int, default=3,
        help='threads [3]')
    args = parser.parse_args()

    # check Flye is installed if graph is not provided
    if not args.input_graph: # check in input is sequences or alignments
        output = subprocess.check_output("{} --version".format(args.flye_exec), shell=True)
    file_extension = os.path.splitext(args.input[0])[1]

    # check minigraph is installed if alignments are not provided
    alignments = []
    if file_extension != ".gaf": # input is sequences, not alignment
        output = subprocess.check_output("{} --version".format(args.minigraph_exec), shell=True)
    else:
        if not args.input_graph: # if alignments are provided, graph must be as well
            raise ValueError("Assembly graph for provided alignments must also be provided.")
        if not args.bp_table: # if alignments are provided, bp_table must be as well
            raise ValueError("Number of bps per sequence must be provided for each alignment.")
        alignments = args.input

    # create output directories
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    if not os.path.isabs(args.output_dir):
        args.output_dir = os.path.normpath(os.path.join(os.getcwd(), args.output_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Flye - create assembly graph
    if not args.input_graph:
        flye_output = os.path.join(args.output_dir, 'flye_output')
        subprocess.check_output("{} {} {} --out-dir {} --threads {} --meta --keep-haplotypes"
                            .format(args.flye_exec, ''.join(["--", args.type]),
                                    " ".join(args.input), flye_output, args.threads),
                                shell=True)
        args.input_graph = os.path.join(flye_output, "assembly_graph.gfa")
        logging.info("Metaflye output complete: %s", flye_output)

    # Minigraph - align sequences to assembly graph
    reads_bps = []
    if file_extension != ".gaf": # input is sequences, not alignment
        for input_seq in args.input:
            file_stem = os.path.splitext(os.path.basename(input_seq))[0]
            if file_extension == ".gz":
                file_stem = os.path.splitext(file_stem)[0]
            alignment_output = os.path.join(args.output_dir, "".join([file_stem, ".gaf"]))
            subprocess.check_output("{} -t{} {} {} > {}"
                                    .format(args.minigraph_exec, args.threads, args.input_graph,
                                            input_seq, alignment_output), shell=True)
            alignments.append(alignment_output)
            reads_bps.append((file_stem, count_bps(input_seq)))
        logging.info("Minigraph alignments complete: %s", args.output_dir)
    # report number of bp per input data set
    reads_bp_df = pd.DataFrame(data=reads_bps)
    reads_bp_out_path = os.path.join(args.output_dir, "bp_counts.tsv")
    reads_bp_df.to_csv(reads_bp_out_path, sep='\t', header=False, index=False)

    # Rhea - evaluate variations in graph alignment coverage
    print(alignments)
