#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""
import os
import re
import csv
import logging
import argparse
import subprocess
from functools import reduce

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
SV_COLUMN_NAMES = ['SV_type', 'neighbor1', 'neighbor2', 'replacement']
NODE_LENGTH_DICT = {}
SEQ_BP_DICT = {}

def convert_edge_dict(in_path):
    """
    Reads in the output matrix of edge weights and converts to a dict of dicts.
        Used for debugging.

    :param in_path: path for .tsv file matrix
    :return: (dict) dict of dicts
    """
    dict_of_dicts = {}
    # Read the TSV file and convert it into a dictionary of dictionaries
    with open(in_path, 'r', encoding='utf-8') as tsvfile:
        tsvreader = csv.DictReader(tsvfile, delimiter='\t')
        for row in tsvreader:
            dict_of_dicts[row['header1']] = {k: v for k, v in row.items() if k != 'header1'}
    return dict_of_dicts

def count_bps(sequence_file):
    """
    Counts the number of bps in a sequence file. Expected format: fastq/a, option .gz compression.

    :str sequence_file: path to input set of sequences
    :return: (int) number of bps in sequence file
    """
    base_file_ext = os.path.splitext(sequence_file.split(".gz")[0])[1]
    regex_command = '/^>/ {{next}}' if base_file_ext[-1] == "a" else 'NR%4==2'
    print_command = 'cat'
    if sequence_file.endswith('.gz'): # check if gzcat or zcat
        print_command = "zcat"
        try:
            subprocess.check_output("{} {}".format(print_command, sequence_file),
                                    shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            print_command = "zcat"
            try:
                subprocess.check_output("{} {}".format(print_command, sequence_file), shell=True)
            except subprocess.CalledProcessError:
                print("Unable to gather bp count from input reads. "
                      "No normalization weights will be applied.")
                return 1
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
        print("Unable to gather bp count from input reads. "
              "No normalization weights will be applied.")
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
                    graph.add_edge(node_a, node_b, count=1)
                else:
                    curr_count = graph[node_a][node_b]["count"]
                    graph.add_edge(node_a, node_b, count=curr_count+1)
    return graph, pd.DataFrame(entries, columns=NODE_DF_HEADERS)

def add_node_coverage_to_graph(graph, df_nodes_coverage):
    """
    Adds a parameter 'log_fold_change_vector' to each node in the networkx graph
        based on the data in df_nodes_coverage.

    :param graph: networkx graph
    :param df_nodes_coverage: dataframe where 'node_id' corresponds to the name
        of the nodes in the graph and [log_fold_change_$t] contains columns with
        log fold change coverage for the node at timepoint $t
    :return: (networkx G) updated networkx graph
    """
    log_fold_change_columns = ['node_id'] + \
                             [col for col in df_nodes_coverage.columns
                              if re.match(r'log_fold_change_\d+$', col)]
    df_nodes_log_fold_change = df_nodes_coverage[log_fold_change_columns].set_index("node_id")
    for node in graph.nodes():
        if node in df_nodes_log_fold_change.index:
            log_fold_change_vector = [float(value) for value
                                     in df_nodes_log_fold_change.loc[node].tolist()]
            graph.nodes[node]['log_fold_change_vector'] = log_fold_change_vector
        else:
            graph.nodes[node]['log_fold_change_vector'] = [0] * df_nodes_log_fold_change.shape[1]
    return graph

def calculate_log_fold_change_edge_coverage(df_nodes):
    """
    Takes in global parameter COVERAGE_EDGE_DICTS and df containing
        a list of all node names under column 'node_id' to produce
        log fold change in edge coverage between subsequent entries
    param df_nodes: df with column 'node_id' containing node names
    :return: (list) list of matricies (dict of dicts) containing
        edge log fold change in coverage.
    """
    # can speed this up with matricies rather than for loops
    lfc_list = []
    for time in range(1, len(COVERAGE_EDGES_DICTS)):
        lfc_matrix = {row: {col: 1 for col in df_nodes['node_id']}
                             for row in df_nodes['node_id']}
        for row_key in COVERAGE_EDGES_DICTS[time]:
            lfc_matrix[row_key] = {}
            for col_key in COVERAGE_EDGES_DICTS[time][row_key]:
                lfc_matrix[row_key][col_key] = np.log2(
                    COVERAGE_EDGES_DICTS[time][row_key][col_key]
                    / COVERAGE_EDGES_DICTS[time-1][row_key][col_key])
        lfc_list.append(lfc_matrix)
    return lfc_list


def label_indel(vectors_stacked, std_thresh=1):
    """
    Takes in the log fold change in edge coverage for a graph triangle
        and outputs locations containing an indel.

    :param vectors_stacked: stacked vectors of log fold change
    :param std_thresh: number of stds away from the median needed to express outlier
    :return: (dict) key is the index location within inner list containing an indel (timepoint)
        value is the index of the vector containing the indel (node)
    """
    # Calculate mean and standard deviation across the three vectors for each index
    medians = np.median(vectors_stacked, axis=0)
    std_devs = np.std(vectors_stacked, axis=0)
    # Identify outliers for each index across the vectors
    outliers = np.abs(vectors_stacked - medians) > std_thresh * std_devs
    # Get indices of outliers
    outlier_indices = np.where(np.sum(outliers, axis=0) == 1)[0]
    outlier_vectors = {index: np.where(outliers[:, index])[0][0] for index in outlier_indices}
    return outlier_vectors

def label_mutations(cycle, vectors_stacked, datas, std_thresh=1):
    """
    Detected mutation in a cycle of length 4 and updated datas dictionary accordingly.

    :param cycle: list of 4 nodes
    :param vectors_stacked: corresponding log fold change coverage vectors stacked
    :param datas: dict of SVs detected
    :param std_thresh: number of stds away from the median needed to express mutation
    :return: (dict) updated datas to contain mutation in this cycle
    """
    medians = np.median(vectors_stacked, axis=0)
    std_devs = np.std(vectors_stacked, axis=0)
    within_threshold = np.abs(vectors_stacked - medians) < std_thresh * std_devs
    mut_13 = np.where(np.all(within_threshold[[0, 2]], axis=0) &
                                ~np.all(within_threshold[[1, 3]], axis=0))[0]
    mut_02 = np.where(np.all(within_threshold[[1, 3]], axis=0) &
                                ~np.all(within_threshold[[0, 2]], axis=0))[0]
    for (timestep_significant, i_mut) in ((mut_13, 1), (mut_02, 0)):
        for timestep in timestep_significant:
            if vectors_stacked[i_mut][timestep] * \
                    vectors_stacked[i_mut+2][timestep] < 0: # one increases and the other decreases
                mut_win, mut_loss = i_mut, i_mut+2
                if vectors_stacked[i_mut][timestep] < vectors_stacked[i_mut+2][timestep]:
                    mut_win, mut_loss = i_mut+2, i_mut
                datas[timestep][cycle[mut_win]].append(
                    ("mutation win", cycle[mut_loss], cycle[i_mut+1], cycle[(i_mut+3)%4]))
                datas[timestep][cycle[mut_loss]].append(
                    ("mutation loss", cycle[mut_win], cycle[i_mut+1], cycle[(i_mut+3)%4]))
    return datas

def detect_structual_variants(graph, df_nodes, node_std_thresh, edge_lfc_thresh):
    """
    Function to detect and return SVs in graph.

    :param graph: networkx graph
    :param df_nodes: df containing node names under column 'node_id'
    :param node_std_thresh: number of stds away from the median needed to express indel or mutation
    :param edge_lfc_thresh: lfc change minimum for change in edge coverage to express duplication
    :return: list of SV data where each entry is a dict of nodes and their SVs for a given timepoint
    """
    datas = [{node: [] for node in df_nodes['node_id']} for _ in range(N_TIMESTEPS)]
    edge_coverage_fold_change = calculate_log_fold_change_edge_coverage(df_nodes)
    cycles = sorted(nx.simple_cycles(graph, length_bound=4))
    for cycle in cycles:
        if len(cycle) == 1:  # look for tandem duplication/deletion
            for timestep in range(N_TIMESTEPS):
                if edge_coverage_fold_change[timestep][cycle[0]][cycle[0]] \
                        > edge_lfc_thresh:
                    datas[timestep][cycle[0]].append(('tandem duplication gain', '', '', ''))
                elif edge_coverage_fold_change[timestep][cycle[0]][cycle[0]] < \
                        (edge_lfc_thresh * -1):
                    datas[timestep][cycle[0]].append(('tandem duplication loss', '', '', ''))
        elif len(cycle) == 3:  # detect indels
            vectors = np.vstack([graph.nodes[node]['log_fold_change_vector'] for node in cycle])
            outlier_vectors = label_indel(vectors, node_std_thresh)
            for timestep, outlier_index in outlier_vectors.items():
                svariant_type = 'deletion' if vectors[outlier_index][timestep] == \
                                              min(vectors[:, timestep]) else 'insertion'
                neighbors = cycle.copy()
                neighbors.pop(outlier_index)
                datas[timestep][cycle[outlier_index]]\
                    .append((svariant_type, neighbors[0], neighbors[1], ''))
        elif len(cycle) == 4: # detect mutation swaps
            vectors = np.vstack([graph.nodes[node]['log_fold_change_vector'] for node in cycle])
            datas = label_mutations(cycle, vectors, datas, node_std_thresh)
    for node_a, node_b, edge_count in graph.edges(data='count', default=1):
        # second possible duplication pattern
        if edge_count >= 2: #look for tandem duplication
            for timestep in range(N_TIMESTEPS):
                if edge_coverage_fold_change[timestep][node_a][node_b] > edge_lfc_thresh:
                    dup_node = node_a
                    if graph.nodes[node_b]['log_fold_change_vector'][timestep] \
                            > graph.nodes[node_a]['log_fold_change_vector'][timestep]:
                        dup_node = node_b
                    if not datas[timestep][dup_node]:
                        datas[timestep][dup_node].append(['tandem duplication gain', '', '', ''])
    return datas

def output_sv_detection_files(datas, df_nodes, out_dir):
    """
    Output a tsv of SV files for each timestep and return a combined collapsed df

    :param datas: list of SV data where each entry is a dict of nodes
        and their SVs for a given timepoint
    :param df_nodes: df containing node names under column 'node_id'
    :param out_dir: path to directory for all output files
    :return: (df) combined collapsed df of SVs for each timepoint
    """
    # add SVs for each timestep into nodes_df
    counter = 0
    merged_sv_dfs = [df_nodes.copy()]
    for data in datas:
        column_names_timestep = ['node_id'] + \
                                ["t{}-{}".format(counter, value) for value in SV_COLUMN_NAMES]
        new_df = pd.DataFrame([(node_id, *tup) for node_id, tuples in data.items()
                               for tup in tuples], columns=column_names_timestep)
        new_merged_df = pd.merge(df_nodes.copy(), new_df, on='node_id', how='left')
        new_df_outpath = os.path.join(out_dir, "structual_variants-t{}.tsv".format(counter))
        new_merged_df.to_csv(new_df_outpath, index=False, sep='\t')
        new_df_collapsed = new_df.groupby('node_id').agg(lambda x: ', '.
                                                         join(filter(None, map(str, set(x)))))
        merged_sv_dfs.append(new_df_collapsed)
        counter += 1
    df_merged = reduce(lambda left, right: pd.merge(left, right, on='node_id', how='outer'),
                       merged_sv_dfs)
    return df_merged

def process_row_coverage(row, timestep):
    """
    Parses row in the .paf alignment file and updates COVERAGE_DICTS
        global parameter.

    :param row: row in the .paf file
    :param timestep: current timestep for row to be processed
    :return: (pd series) containing 'node tuples' entry of each coveraged node
        and respective bps covered
    """
    # gather list of all nodes in alignment
    target_name_parts = [item for item in re.split('<|>| ', str(row['target name'])) if item]
    if len(target_name_parts) == 1: # if alignment is only to one node
        node_align_length = row['target end'] - row['target start']
        part_trim = target_name_parts[0].split("_")[-1]
        COVERAGE_DICTS[timestep][part_trim] += node_align_length

    # measure number of bp in each node included in alignment
    else:
        for j, part in enumerate(target_name_parts):
            part_trim = part.split("_")[-1]
            node_length = NODE_LENGTH_DICT[part]
            if j == 0:
                node_align_length = node_length - row['target start']
            elif j == len(target_name_parts) - 1:
                node_align_length = node_length - (row['target length'] - row['target end'])
            else:
                node_align_length = node_length
            COVERAGE_DICTS[timestep][part_trim] = COVERAGE_DICTS[timestep][part_trim] + \
                                                  node_align_length
            if not j == 0: # add to edge coverage
                COVERAGE_EDGES_DICTS[timestep][part_trim][prev_part_trim] += 1
                if part_trim != prev_part_trim:
                    COVERAGE_EDGES_DICTS[timestep][prev_part_trim][part_trim] += 1
            prev_part_trim = part_trim


def create_coverage_df(alignments_pafs, df_nodes):
    """
    Add to df_nodes with the coverage of each node based on the supplied alignments_pafs.

    :param alignments_gafs: list of paths of alignment .gaf in order of the series
    :return: (df) df nodes and (list) coverage cols of column names containing new coverage values
    """
    # for each alignment, measure node coverage
    coverage_cols = []
    timestep = 0

    for alignment_paf in alignments_pafs:
        align_stem = os.path.splitext(os.path.basename(alignment_paf))[0]
        logging.info("Calculating graph coverage for: %s", align_stem)
        paf_df = pd.read_csv(alignment_paf, names=PAF_HEADERS, sep='\t')
        # Apply the function to each row to update node coverage lengths
        paf_df.apply(process_row_coverage, axis=1, timestep=timestep)

        # add updated coverage to df_nodes, normalizing for the length of the node
        node_coverage_len_df = pd.DataFrame(list(COVERAGE_DICTS[timestep].items()),
                                            columns=['node_id', 'coverage_bps'])
        df_nodes = pd.merge(df_nodes, node_coverage_len_df, on='node_id', how='left')
        df_nodes[align_stem] = df_nodes['coverage_bps'] / df_nodes['node_length']
        df_nodes = df_nodes.drop(columns=['coverage_bps'])
        coverage_cols.append(align_stem)
        timestep += 1
    return df_nodes, coverage_cols

def create_coverage_normalized_df(df_nodes_coverage, coverage_cols):
    """
    Adjust each of the coverage_cols in df_nodes_coverage to be nomalized based on the number of
        base pairs in the aligned sample.

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
        group_max with seaborn color palette.

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
    Add colors values to coverage_cols in df_nodes_coverages.

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
                                         normalized=False, lfc=False):
    """
    Add colors values to diff_coverage_cols in df_nodes_coverages to show change
        between elements in series data.

    :param df_nodes_coverage: dataframe with coverage values
    :param diff_coverage_cols: columns containing coverage values difference
        between consecutive series elements
    :param normalized: boolean of if these are normalized values.
        Only used for output string.
    :param lfc: boolean of if these are log fold change (vs linear) difference values.
        Only used for output string.
    :return: (df) updated with color hex values to correspond to coverage diff values.
    """
    coverage_file = " normalized" if normalized else ""
    lfc_or_linear = "log fold change" if lfc else "linear change"

    df_nodes_coverage = df_nodes_coverage.fillna(0)
    max_diff = max(df_nodes_coverage[diff_coverage_cols].quantile(.95))
    min_diff = max(df_nodes_coverage[diff_coverage_cols].quantile(.05))
    min_max = max(max_diff, abs(min_diff))
    logging.info("Min/Max val for coverage%s %s difference color spectrum: %s",
                 coverage_file, lfc_or_linear, min_max)
    for col in diff_coverage_cols:
        df_nodes_coverage[f"{col}_colour"] = df_nodes_coverage[col].apply(map_to_heatmap_color,
                                                                      group_max=min_max,
                                                                      group_min=(min_max * -1),
                                                                      palette="coolwarm")
    return df_nodes_coverage

def calculate_diff_to_nodes_df(df_nodes_coverage, coverage_cols):
    """
    Add 2 columns for each consecutive col in coverage_cols: 1 for
        linear difference and 1 for lfc difference.

    :param df_nodes_coverage: dataframe with coverage values
    :param coverage_cols: list of column names containing coverage information
    :return: (df) updated df with consecutive difference columns, (list) diff_cols
        of columns containing linear difference, (list) diff_lfc_cols of columns
        containing lfc difference
    """
    linear_change_cols, log_fold_change_cols = [], []
    # set all coverage below 1 to 1 for calculating difference
    df_coverage_adjusted = pd.DataFrame()
    df_coverage_adjusted[coverage_cols] = \
        df_nodes_coverage[coverage_cols].applymap(lambda x: 1 if x < 1 else x)
    for counter in range(len(coverage_cols)-1):
        pair = (coverage_cols[counter], coverage_cols[counter+1])
        new_linear_change_col = f"linear_change_{counter}"
        df_nodes_coverage[new_linear_change_col] = df_coverage_adjusted[f"{pair[1]}"] \
                                                   - df_coverage_adjusted[f"{pair[0]}"]
        linear_change_cols.append(new_linear_change_col)
        new_log_fold_change_col = f"log_fold_change_{counter}"
        df_nodes_coverage[new_log_fold_change_col] = np.log2(df_coverage_adjusted[pair[1]]
                                                             / df_coverage_adjusted[pair[0]])
        log_fold_change_cols.append(new_log_fold_change_col)
    return df_nodes_coverage, linear_change_cols, log_fold_change_cols

def create_both_nodes_coverage_dfs(alignments_list_ordered, df_nodes):
    """
    Generate both df of nodes coverage and a normalized version. Contains coverage
        for each alignment in alighments_list_ordered with corresponding hex color values.
        Contains value for differnce between consectuive alignments (linear and lfc)
        with corresponding hex values.

    :param alignments_list_ordered: list of paths to alignments
    :param df_nodes: dataframe containing all nodes in graph
    :return: (list) of two dataframes: raw and normalized
    """
    df_nodes_coverage, alignment_stems_list = create_coverage_df(alignments_list_ordered, df_nodes)
    df_nodes_coverage_norm = create_coverage_normalized_df(df_nodes_coverage, alignment_stems_list)
    output_dfs = []
    for (df_node, normalized_flag) in [(df_nodes_coverage, False), (df_nodes_coverage_norm, True)]:
        df_node = add_coverage_colors_to_nodes_df(df_node, alignment_stems_list,
                                                  normalized=normalized_flag)
        df_node, cols_change, cols_change_lfc = calculate_diff_to_nodes_df(
            df_node, alignment_stems_list)
        df_node = add_coverage_diff_colors_to_nodes_df(df_node, cols_change,
                                                       normalized=normalized_flag)
        df_node = add_coverage_diff_colors_to_nodes_df(df_node, cols_change_lfc,
                                                  normalized=normalized_flag, lfc=True)
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
    SEQ_BP_DICT = dict(zip(reads_bp_df[0], reads_bp_df[1]))
    for in_file in alignments:
        in_file_stem = os.path.splitext(os.path.basename(in_file))[0]
        if in_file_stem not in SEQ_BP_DICT:
            raise AssertionError("{} is not included in supplied bp_table as expected"
                                 .format(in_file_stem))

    # Start Rhea - read in graph & init variables
    networkx_graph, nodes_df = read_graph(args.input_graph)
    NODE_LENGTH_DICT = dict(zip(nodes_df['node'], nodes_df['node_length']))
    N_INPUTS = len(args.input)
    COVERAGE_DICTS = [{key: 0 for key in nodes_df['node_id']} for _ in range(N_INPUTS)]
    COVERAGE_EDGES_DICTS = [{row: {col: 1 for col in nodes_df['node_id']}
                             for row in nodes_df['node_id']} for _ in range(N_INPUTS)]
    N_TIMESTEPS = N_INPUTS - 1

    # calculate edge coverage and log fold change
    nodes_df_coverage, nodes_df_coverage_norm = \
        create_both_nodes_coverage_dfs(alignments, nodes_df)

    # output edge_coverage
    for i, in_file in enumerate(alignments):
        file_stem = os.path.splitext(os.path.basename(in_file))[0]
        df_edge_coverage = pd.DataFrame.from_dict(COVERAGE_EDGES_DICTS[i], orient='index')
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
    variants_data = detect_structual_variants(networkx_graph, nodes_df,
                                            args.node_std, args.edge_lfc_thresh)
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
