#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""
import re
import os
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns



__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'Jan 2024'

NODE_DF_HEADERS = ['node', 'node_id', 'node_length']
PAF_HEADERS = ['query', 'query length', 'query start', 'query end', 'relative', 'target name',
               'target length', 'target start', 'target end', 'n matches', 'alignment length',
               'quality', 'tp', 'cm', 's1', 's2', 'dv']
NODE_LENGTH_DICT = {}
SEQ_BP_DICT = {}


def count_bps(sequence_file):
    """
    Counts the number of bps in a sequence file. Expected format: fastq/a, option .gz compression.

    :str sequence_file: path to input set of sequences
    :return: (int) number of bps in sequence file
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


def read_graph(graph_path):
    """
    Read in .gfa Metaflye graph to a networkx graph and
        create a df of all graph nodes and their length in bps

    :str graph_path: path to graph.gfa file
    :return: (df) of graph nodes with their length in bps
    """

    logging.info("Reading in graph: %s", graph_path)
    graph = nx.Graph()
    entries = []
    with open(graph_path, "r", encoding='UTF-8') as file:
        for line in file:
            if line.startswith("S"): # add a node
                split_line = line[:-1].split("\t")
                new_node = split_line[1]
                new_node_id = split_line[1].split("_")[-1]
                node_length = len(split_line[2])
                entries.append((new_node, new_node_id, node_length))
            if line.startswith("L"): # add an edge
                split_line = line.rstrip("\n").split("\t")
                node_a = split_line[1].split("_")[-1]
                node_b = split_line[3].split("_")[-1]
                if not graph.has_edge(node_a, node_b):
                    graph.add_edge(node_a, node_b)
    return graph, pd.DataFrame(entries, columns=NODE_DF_HEADERS)

def add_node_coverage_to_graph(graph, df_nodes_coverage):
    """

    :param graph:
    :param df_nodes_coverage:
    :return:
    """
    percent_change_columns = ['node_id'] + \
                             [col for col in df_nodes_coverage.columns if re.match(r'change_percent_\d+$', col)]
    df_nodes_percent_change = df_nodes_coverage[percent_change_columns].set_index("node_id")
    for node in graph.nodes():
        if node in df_nodes_percent_change.index:
            percent_change_vector = [float(value) for value in df_nodes_percent_change.loc[node].tolist()]
            graph.nodes[node]['percent_change_vector'] = percent_change_vector
        else:
            graph.nodes[node]['percent_change_vector'] = [0] * df_nodes_percent_change.shape[1]
    return graph

def label_indel(vectors):
    num_vectors = len(vectors)
    vector_length = len(vectors[0])
    # Combine the vectors into a single 2D array
    combined = np.vstack(vectors)
    # Calculate mean and standard deviation across the three vectors for each index
    median = np.median(combined, axis=0)
    std_dev = np.std(combined, axis=0)
    # Define a threshold for outliers (e.g., values more than 3 standard deviations away from the mean)
    threshold = 1  # Adjust as needed based on your data and requirements
    # Identify outliers for each index across the vectors
    outliers = np.abs(combined - median) > threshold * std_dev
    # Get indices of outliers
    outlier_indices = np.where(np.sum(outliers, axis=0) == 1)[0]
    outlier_vectors = {index: np.where(outliers[:, index])[0][0] for index in outlier_indices}
    return outlier_vectors


def detect_structual_variants(graph, df_nodes):
    """

    :param graph:
    :param df_nodes:
    :return:
    """
    df_SVs = df_nodes.copy()
    datas = [[]] * N_TIMESTEPS
    cliques = list(nx.enumerate_all_cliques(graph))
    for clique in cliques:
        if len(clique) == 3: # detect indels
            vectors = [graph.nodes[node]['percent_change_vector'] for node in clique]
            outlier_vectors = label_indel(vectors)
            for timestep, outlier_index in outlier_vectors.items():
                SV_type = 'deletion' if vectors[outlier_index] == min(vectors) else 'insertion'
                neighbors = clique.copy()
                neighbors.pop(outlier_index)
                datas[timestep].append([clique[outlier_index], SV_type, neighbors[0], neighbors[1], ""])

    # add SVs for each timestep into nodes_df
    columns_names = ['SV_type', 'neighbor1', 'neighbor2', 'replacement']
    counter = 1
    for data in datas:
        column_names_timestep = ['node_id'] + ["t{}-{}".format(counter, value) for value in columns_names]
        additional_df = pd.DataFrame(data=data, columns=column_names_timestep)
        df_SVs = pd.merge(df_SVs, additional_df, on='node_id', how='outer')
        counter = counter + 1
    return df_SVs

def process_row_coverage(row):
    """
    helper function to parse each row in the .paf alignment file

    :param row: row in the .paf file
    :return: (pd series) containing 'node tuples' entry of each coveraged node
        and respective bps covered
    """
    # gather list of all nodes in alignment
    target_name_parts = [item for item in re.split('<|>| ', row['target name']) if item]
    alignment_node_length = []

    if len(target_name_parts) == 1: # if alignment is only to one node
        node_align_length = row['target end'] - row['target start']
        alignment_node_length.append([target_name_parts[0], node_align_length])

    # measure number of bp in each node included in alignment
    else:
        for i, part in enumerate(target_name_parts):
            node_length = NODE_LENGTH_DICT[part]
            if i == 0:
                node_align_length = node_length - row['target start']
            elif i == len(target_name_parts) - 1:
                node_align_length = node_length - (row['target length'] - row['target end'])
            else:
                node_align_length = node_length
            alignment_node_length.append([part, node_align_length])

        # return series containing
    return pd.Series({
        'query': row['query'],
        'query start': row['query start'],
        'query end': row['query end'],
        'relative': row['relative'],
        'target name': row['target name'],
        'target length': row['target length'],
        'target start': row['target start'],
        'target end': row['target end'],
        'alignment length': row['alignment length'],
        'quality': row['quality'],
        'node tuples': alignment_node_length})

def create_coverage_df(alignments_pafs, df_nodes):
    """
    add to df_nodes with the coverage of each node based on the supplied alignments_pafs

    :param alignments_gafs: list of paths of alignment .gaf in order of the series
    :return: (df) df nodes and (list) coverage cols of column names containing new coverage values
    """
    # for each alignment, measure node coverage
    coverage_cols = []
    for alignment in alignments_pafs:
        alignment_stem = os.path.splitext(os.path.basename(alignment))[0]
        paf_df = pd.read_csv(alignment, names=PAF_HEADERS, sep='\t')

        # Apply the function to each row and concatenate the results
        temp_df = paf_df.apply(process_row_coverage, axis=1)

        # Explode the target_name and alignment_node_length lists into separate rows
        temp_df = temp_df.explode('node tuples')
        df_coverage = temp_df['node tuples'].\
            apply(lambda x: pd.Series(x, index=['node', 'node alignment length']))
        df_coverage = df_coverage.reset_index(drop=True)
        df_coverage = df_coverage.groupby('node').sum().reset_index()
        df_coverage['node length'] = df_coverage['node'].map(NODE_LENGTH_DICT)
        df_coverage['coverage'] = df_coverage['node alignment length'] / df_coverage['node length']

        # add coverage to df_nodes
        node_coverage_dict = dict(zip(df_coverage['node'], df_coverage['coverage']))
        df_nodes[alignment_stem] = df_nodes['node'].map(node_coverage_dict)
        coverage_cols.append(alignment_stem)
    return df_nodes, coverage_cols

def create_coverage_normalized_df(df_nodes_coverage, coverage_cols):
    """
    adjust each of the coverage_cols in df_nodes_coverage to be nomalized based on the number of
    base pairs in the aligned sample

    :param df_nodes_coverage: df containing coverage of each node
    :param coverage_cols: column headers for the columns containing series coverage information
    :return: (df) of normalized coverage values
    """
    # add coverage for normalized by total bps per sequence
    df_nodes_coverage_norm = df_nodes_coverage.copy()[NODE_DF_HEADERS]
    median = np.median(list(SEQ_BP_DICT.values()))
    weights = {key: median / value for key, value in SEQ_BP_DICT.items()}
    for col in coverage_cols:
        df_nodes_coverage_norm[col] = df_nodes_coverage[col] * weights[col]
    return df_nodes_coverage_norm

def map_to_heatmap_color(value, group_max, group_min=0, palette="rocket_r"):
    """
    Generates color hex value for value on a heatmap scale from group_min to
        group_max with seaborn color palette

    :param value: float value to calculate heatmap color of
    :param group_max: float max for scale
    :param group_min: float min for scale
    :param palette: sns color palatte
    :return: (str) hex color value
    """
    # seaborn's color_palette is used to create a heatmap-like color spectrum
    cmap = sns.color_palette(palette, as_cmap=True)
    normalized_value = (value - group_min) / (group_max - group_min)
    rgb_color = cmap(normalized_value)[:3]  # Extract RGB components
    hex_color = "#{:02X}{:02X}{:02X}".format(int(255 * rgb_color[0]),
                                             int(255 * rgb_color[1]),
                                             int(255 * rgb_color[2]))
    return hex_color

def add_coverage_colors_to_nodes_df(df_nodes_coverage, coverage_cols, normalized=False):
    """
    Add colors values to coverage_cols in df_nodes_coverages

    :param df_nodes_coverage: dataframe with coverage values
    :param coverage_cols: columns containing coverage values
    :param normalized: boolean of if these are normalized values. Only used for output string.
    :return: updated df with color hex values to correspond to coverage values.
    """
    coverage_file = ""
    if normalized:
        coverage_file = " normalized"
    df_nodes_coverage = df_nodes_coverage.fillna(0)
    max_cov = max(df_nodes_coverage[coverage_cols].quantile(.98))
    logging.info("Maximum val for coverage%s color spectrum: %s",
                 coverage_file, max_cov)
    for col in coverage_cols:
        df_nodes_coverage[f"{col}_colour"] = df_nodes_coverage[col].apply(map_to_heatmap_color,
                                                                      group_max=max_cov)
    return df_nodes_coverage

def add_coverage_diff_colors_to_nodes_df(df_nodes_coverage, diff_coverage_cols,
                                         normalized=False, percentage=False):
    """
    Add colors values to diff_coverage_cols in df_nodes_coverages to show change
        between elements in series data.

    :param df_nodes_coverage: dataframe with coverage values
    :param diff_coverage_cols: columns containing coverage values difference
        between consecutive series elements
    :param normalized: boolean of if these are normalized values.
        Only used for output string.
    :param percentage: boolean of if these are percentage (vs linear) difference values.
        Only used for output string.
    :return: updated df with color hex values to correspond to coverage diff values.
    """
    coverage_file = " normalized" if normalized else ""
    percent_or_diff = "percent" if percentage else "linear"

    df_nodes_coverage = df_nodes_coverage.fillna(0)
    max_diff = max(df_nodes_coverage[diff_coverage_cols].quantile(.95))
    min_diff = max(df_nodes_coverage[diff_coverage_cols].quantile(.05))
    min_max = max(max_diff, abs(min_diff))
    logging.info("Min/Max val for coverage%s %s difference color spectrum: %s",
                 coverage_file, percent_or_diff, min_max)
    for col in diff_coverage_cols:
        df_nodes_coverage[f"{col}_colour"] = df_nodes_coverage[col].apply(map_to_heatmap_color,
                                                                      group_max=min_max,
                                                                      group_min=(min_max * -1),
                                                                      palette="coolwarm")
    return df_nodes_coverage



def calculate_diff_to_nodes_df(df_nodes_coverage, coverage_cols):
    """
    Add 2 columns for each consecutive col in coverage_cols: 1 for
        linear difference and 1 for percentage difference

    :param df_nodes_coverage: dataframe with coverage values
    :return: (df) updated df with consecutive difference columns, (list) diff_cols
        of columns containing linear difference, (list) diff_percent_cols of columns
        containing percent difference
    """
    diff_cols = []
    diff_percent_cols = []
    for counter in range(len(coverage_cols)-1):
        pair = (coverage_cols[counter], coverage_cols[counter+1])
        new_col_diff = f"change_{counter}"
        df_nodes_coverage[new_col_diff] = df_nodes_coverage[f"{pair[1]}"] \
                                          - df_nodes_coverage[f"{pair[0]}"]
        diff_cols.append(new_col_diff)
        new_col_diff_percent = f"change_percent_{counter}"
        df_nodes_coverage[new_col_diff_percent] = df_nodes_coverage[new_col_diff]/\
                                                  df_nodes_coverage[f"{pair[0]}"].replace(0, 0.001)
        diff_percent_cols.append(new_col_diff_percent)
    return df_nodes_coverage, diff_cols, diff_percent_cols

def create_both_nodes_coverage_dfs(alignments_list_ordered, df_nodes):
    """
    Generate both df of nodes coverage and a normalized version. Contains coverage
        for each alignment in alighments_list_ordered with corresponding hex color values.
        Contains value for differnce between consectuive alignments (linear and percentage)
        with corresponding hex values.

    :param alignments_list_ordered: list of paths to alignments
    :param df_nodes: dataframe containing all nodes in graph
    :return: list of two dataframes: raw and normalized
    """
    df_nodes_coverage, alignment_stems_list = create_coverage_df(alignments_list_ordered, df_nodes)
    df_nodes_coverage_norm = create_coverage_normalized_df(df_nodes_coverage, alignment_stems_list)
    output_dfs = []
    for (df_node, normalized_flag) in [(df_nodes_coverage, False), (df_nodes_coverage_norm, True)]:
        df_node = add_coverage_colors_to_nodes_df(df_node, alignment_stems_list,
                                                  normalized=normalized_flag)
        df_node, cols_change, cols_change_percent = calculate_diff_to_nodes_df(
            df_node, alignment_stems_list)
        df_node = add_coverage_diff_colors_to_nodes_df(df_node, cols_change,
                                                       normalized=normalized_flag)
        df_node = add_coverage_diff_colors_to_nodes_df(df_node, cols_change_percent,
                                                  normalized=normalized_flag, percentage=True)
        output_dfs.append(df_node)
    return output_dfs[0], output_dfs[1]


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
        '--threads', '-t', type=int, default=3,
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
        flye_output = os.path.join(args.output_dir, 'flye_output')
        subprocess.check_output("{} {} {} --out-dir {} --threads {} --meta --keep-haplotypes"
                            .format(args.flye_exec, ''.join(["--", args.type]),
                                    " ".join(args.input), flye_output, args.threads),
                                shell=True)
        args.input_graph = os.path.join(flye_output, "assembly_graph.gfa")
        logging.info("Metaflye output complete: %s", flye_output)
    else:
        if not os.path.exists(args.input_graph):
            raise FileNotFoundError(f"The file at {args.input_graph} does not exist.")
        logging.info("Reading in Metaflye graph from: %s", args.input_graph)

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
            logging.info("Reading in alignment paf: %s", alignment_gaf)


    # Rhea - evaluate variations in graph alignment coverage
    networkx_graph, nodes_df = read_graph(args.input_graph)
    NODE_LENGTH_DICT = dict(zip(nodes_df['node'], nodes_df['node_length']))
    SEQ_BP_DICT = dict(zip(reads_bp_df[0], reads_bp_df[1]))
    N_TIMESTEPS = len(args.input) - 1
    nodes_df_coverage, nodes_df_coverage_norm = \
        create_both_nodes_coverage_dfs(alignments, nodes_df)

    #nodes_df_coverage = pd.read_csv("/Users/kcurry/hotspring-metagenomic/RHEA/m0_rhea/node_coverage.csv", dtype=str)
    #nodes_df_coverage_norm = pd.read_csv("/Users/kcurry/hotspring-metagenomic/RHEA/m0_rhea/node_coverage_norm.csv", dtype=str)

    networkx_graph = add_node_coverage_to_graph(networkx_graph, nodes_df_coverage)
    variants_df = detect_structual_variants(networkx_graph, nodes_df)

    # export coverage dataframes
    variants_df_outpath = os.path.join(args.output_dir, "structual_variants.tsv")
    coverage_df_outpath = os.path.join(args.output_dir, "node_coverage.csv")
    coverage_df_norm_outpath = os.path.join(args.output_dir, "node_coverage_norm.csv")
    variants_df.to_csv(variants_df_outpath, index=False, sep='\t')
    nodes_df_coverage.to_csv(coverage_df_outpath, index=False)
    nodes_df_coverage_norm.to_csv(coverage_df_norm_outpath, index=False)

    print("done")
