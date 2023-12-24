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

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'Jan 2024'

PAF_headers = ['query', 'query length', 'query start', 'query end', 'relative', 'target name',
               'target length', 'target start', 'target end', 'n matches', 'alignment length',
               'quality', 'tp', 'cm', 's1', 's2', 'dv']



def weight_graph():
    """
    Parses .gfa file from MetaFlye output into a networkx graph. Nodes are weighted by length of
    sequence. edges are weighted by coverage.

    :str graph_path: path to .gfa MetaFlye output
    :return: Weighted Networkx graph (G), dictionary of edge id to SeqRecord sequence (sequences)
    """
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    parser.add_argument(
        "input", type=str, nargs='+',
        help="path to metagenome sequences files to compare, in order")
    parser.add_argument(
        "--input_graph", type=str,
        help="path to .gfa assembly graph by MetaFlye")
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
        try: # check flye is installed
            output = subprocess.check_output("{} --version".format(args.flye_exec), shell=True)
        except subprocess.CalledProcessError as e:
            print("Error: Flye not installed or linked: ", e)
    file_extension = os.path.splitext(args.input[0])[1]

    # check minigraph is installed if alignments are not provided
    alignments = []
    if file_extension != ".gaf": # input is sequences, not alignment
        try: # check minigraph is installed
            output = subprocess.check_output("{} --version".format(args.minigraph_exec), shell=True)
        except subprocess.CalledProcessError as e:
            print("Error: Minigraph not installed or linked: ", e)
    else:
        if not args.input_graph: # if alignments are provided, graph must be as well
            raise ValueError("Assembly graph for provided alignments must also be provided.")
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
    if file_extension != ".gaf": # input is sequences, not alignment
        for input_seq in args.input:
            file_stem = os.path.splitext(os.path.basename(input_seq))[0]
            alignment_output = os.path.join(args.output_dir, "".join([file_stem, ".gaf"]))
            subprocess.check_output("{} -t{} {} {} > {}"
                                    .format(args.minigraph_exec, args.threads, args.input_graph,
                                            input_seq, alignment_output), shell=True)
            alignments.append(alignment_output)
        logging.info("Minigraph alignments complete: %s", args.output_dir)

    # Rhea - evaluate variations in graph alignment coverage
    print(alignments)
