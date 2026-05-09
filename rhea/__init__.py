"""
rhea: reference-free heterogeneity and evolution in assembly graphs

Rhea is a software used to detect structural variants (SVs) between steps 
in long-read metagenomic series data.
"""

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'Jan 2024'

import os
import csv
import logging
import argparse
import subprocess

import pandas as pd

# Import all core functions
from .core import (
    NODE_DF_HEADERS, PAF_HEADERS, SV_COLUMN_NAMES, # global variables
    convert_edge_dict, count_bps, read_graph, process_row_coverage, 
    create_coverage_df, create_coverage_normalized_df, map_to_heatmap_color,
    add_coverage_colors_to_nodes_df, calculate_diff_to_nodes_df,
    create_both_nodes_coverage_dfs, add_coverage_diff_colors_to_nodes_df,
    calculate_log_fold_change_edge_coverage, add_node_coverage_to_graph,
    detect_structural_variants, output_sv_detection_files
)

def main():
    
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
        "--node-std", type=float, default=1.0,
        help="number of stds away from the median for log fold change")
    parser.add_argument(
        "--edge-lfc-thresh", type=float, default=1.0,
        help="min value of log fold change for change in edge coverage significance")
    parser.add_argument(
        "--flye-exec", type=str, default='flye',
        help="path to flye executable")
    parser.add_argument(
        "--minigraph-exec", type=str, default='minigraph',
        help="path to minigraph executable")
    parser.add_argument(
        '--output-dir', '-o', type=str, default=os.path.join(os.getcwd(), "rhea_results"),
        help='output directory name [./rhea_results]')
    parser.add_argument(
        '--collapse', action="store_true",
        help='Does not use --keep-haplotypes strain variation for metaFlye')
    parser.add_argument(
        '--raw-diff', action="store_true",
        help='Use raw coverage difference rather than normalized')
    parser.add_argument(
        '--threads', '-t', type=int, default=3,
        help='threads [3]')
    args = parser.parse_args()

    # check more than one graph input is provied:
    if len(args.input) < 2:
        raise ValueError("A minimum of 2 input files are required.")

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
    if not os.path.isabs(args.output_dir):
        args.output_dir = os.path.normpath(os.path.join(os.getcwd(), args.output_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    log_outpath = os.path.join(args.output_dir, 'rhea.log')
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO, handlers=[logging.FileHandler(log_outpath),
                                   logging.StreamHandler()])

    # Flye - create assembly graph
    if not args.input_graph:
        graph_detail = "" if args.collapse else "--keep-haplotype"
        flye_output = os.path.join(args.output_dir, 'metaflye')
        subprocess.check_output("{} {} {} --out-dir {} --threads {} --meta {}"
                            .format(args.flye_exec, ''.join(["--", args.type]),
                                    " ".join(args.input), flye_output, args.threads,
                                    graph_detail),
                                shell=True)
        args.input_graph = os.path.join(flye_output, "assembly_graph.gfa")
        logging.info("Metaflye output complete: %s", flye_output)
    else:
        if not os.path.exists(args.input_graph):
            raise FileNotFoundError(f"The file at {args.input_graph} does not exist.")

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
    else: # if alignments, then read in bp count table
        if not os.path.exists(args.bp_table):
            raise FileNotFoundError(f"The file at {args.bp_table} does not exist.")
        logging.info("Reading in bp_table: %s", args.bp_table)
        reads_bp_df = pd.read_csv(args.bp_table, sep='\t', header=None)
        for alignment_gaf in args.input:
            if not os.path.exists(alignment_gaf):
                raise FileNotFoundError(f"The file at {alignment_gaf} does not exist.")
            if os.path.splitext(alignment_gaf)[1] != ".gaf":
                raise FileNotFoundError(f"The file at {alignment_gaf} is not in format: gaf.")

    # check all sequence files are included in reads_bp_df:
    seq_bp_dict = dict(zip(reads_bp_df[0], reads_bp_df[1]))
    for in_file in alignments:
        in_file_stem = os.path.splitext(os.path.basename(in_file))[0]
        if in_file_stem not in seq_bp_dict:
            raise AssertionError("{} is not included in supplied bp_table as expected"
                                 .format(in_file_stem))

    # Start Rhea - read in graph & init variables
    networkx_graph, nodes_df = read_graph(args.input_graph)
    NODE_LENGTH_DICT = dict(zip(nodes_df['node'], nodes_df['node_length']))
    N_INPUTS = len(args.input)
    N_TIMESTEPS = N_INPUTS - 1

    # calculate edge coverage and log fold change
    nodes_df_coverage, nodes_df_coverage_norm, coverage_dicts, coverage_edges_dicts = \
        create_both_nodes_coverage_dfs(alignments, nodes_df, NODE_LENGTH_DICT, seq_bp_dict)

    # output edge_coverage
    for i, in_file in enumerate(alignments):
        file_stem = os.path.splitext(os.path.basename(in_file))[0]
        df_edge_coverage = pd.DataFrame.from_dict(coverage_edges_dicts[i], orient='index')
        df_edge_coverage_outpath = os.path.join(args.output_dir,
                                                "edge_coverage-{}.tsv".format(file_stem))
        df_edge_coverage.to_csv(df_edge_coverage_outpath, sep='\t')

    # output node coverage
    coverage_df_outpath = os.path.join(args.output_dir, "node_coverage.csv")
    coverage_df_norm_outpath = os.path.join(args.output_dir, "node_coverage_norm.csv")
    nodes_df_coverage.to_csv(coverage_df_outpath, index=False)
    nodes_df_coverage_norm.to_csv(coverage_df_norm_outpath, index=False)

    # detect structural variants
    if args.raw_diff:
        networkx_graph = add_node_coverage_to_graph(networkx_graph, nodes_df_coverage)
    else:
        networkx_graph = add_node_coverage_to_graph(networkx_graph, nodes_df_coverage_norm)
    variants_data = detect_structural_variants(networkx_graph, nodes_df,
                                            args.node_std, args.edge_lfc_thresh,
                                            N_TIMESTEPS, coverage_edges_dicts)
    variants_df = output_sv_detection_files(variants_data, nodes_df, args.output_dir)
    if args.raw_diff:
        complete_df = pd.merge(variants_df, nodes_df_coverage,
                               on=['node', 'node_id', 'node_length'], how='left')
    else:
        complete_df = pd.merge(variants_df, nodes_df_coverage_norm,
                               on=['node', 'node_id', 'node_length'], how='left')
    complete_df_outpath = os.path.join(args.output_dir, "Bandage_metadata.csv")
    complete_df.to_csv(complete_df_outpath, index=False)
    logging.info("Rhea complete: %s", args.output_dir)


if __name__ == "__main__":
    main()

