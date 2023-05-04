#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""

import os
import math
import logging
import argparse
import subprocess
from collections import defaultdict

import pandas as pd
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import SpectralClustering

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'May 2023'


def weight_graph(graph_path):
    """
    Parses .gfa file from MetaFlye output into a networkx graph. Nodes are weighted by length of
    sequence. edges are weighted by coverage.

    :str graph_path: path to .gfa MetaFlye output
    :return: Weighted Networkx graph (G), dictionary of edge id to SeqRecord sequence (sequences)
    """
    logging.info("Reading in graph: %s", graph_path)
    graph = nx.Graph()
    sequences = {}
    with open(graph_path, "r", encoding='UTF-8') as file:

        for line in file:
            if line.startswith("L"):
                split_line = line.rstrip("\n").split("\t")
                node_a = split_line[1].split("_")[-1]
                node_b = split_line[3].split("_")[-1]
                read_coverage = [cov.split(":")[-1] for cov in split_line if cov.startswith("RC")]
                weight = 0
                if read_coverage:
                    weight = int(read_coverage[0])
                # add weight to existing edge weight
                if not graph.has_edge(node_a, node_b):
                    graph.add_edge(node_a, node_b, weight=weight)
                else:
                    weight = graph[node_a][node_b]["weight"] + weight
                    graph.add_edge(node_a, node_b, weight=weight)

            if line.startswith("S"):
                split_line = line[:-1].split("\t")
                new_node = split_line[1].split("_")[-1]
                node_length = len(split_line[2])
                graph.add_node(new_node, weight=node_length)
                sequences[new_node] = SeqRecord(Seq(split_line[2]),
                                                id=new_node,
                                                description="edge_{}".format(new_node))
    return graph, sequences


def cluster_graph(graph, sequences, genome_avg_len, output_path):
    """
    Performs spectral clustering on weighted graph, where number of clusters decided such that
    each cluster is the length of one genome.

    :networkx.G graph: weighted assembly graph. Nodes weighted by length of sequence.
            Edges weighted by coverage between nodes.
    :dict{str:SeqRecord} sequences: dictionary of edge id to SeqRecord sequence
    :int genome_avg_len: expected average genome length for genomes in metagenome sample
    :str output_path: path to output clustered fasta files
    :return: Pandas df mapping graph nodes to assigned clusters
    """
    logging.info("Clustering graphs, depending on the complexity of your graph, "
                 "this may take several minutes; output to: %s", output_path)
    nodes_list, labels_list = [], []
    total_bins = 0
    for connected_component in nx.connected_components(graph):
        adj_matrix = nx.to_numpy_array(graph, nodelist=list(connected_component))
        length_sum = sum([graph.nodes("weight")[node] for node in connected_component])
        n_clusters = math.ceil(length_sum / genome_avg_len)
        if n_clusters == 1:
            next_labels = [total_bins] * len(connected_component)
            cluster_seqs = [sequences.pop(node) for node in connected_component]
            SeqIO.write(cluster_seqs, os.path.join(output_path, "{}.fna".format(total_bins)),
                        "fasta")
        else:
            spec_clust = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
            spec_clust.fit(adj_matrix)
            next_labels = [i + total_bins for i in list(spec_clust.labels_)]
            clusters_dict = defaultdict(list)
            for label, component_nodes in zip(next_labels, list(connected_component)):
                clusters_dict[label].append(component_nodes)
            for k, cluster_nodes in clusters_dict.items():
                cluster_seqs = [sequences.pop(node) for node in cluster_nodes]
                SeqIO.write(cluster_seqs, os.path.join(output_path, "{}.fna".format(k)), "fasta")
        nodes_list = nodes_list + list(connected_component)
        labels_list = labels_list + next_labels
        total_bins = total_bins + n_clusters
    dataframe = pd.DataFrame({'node': nodes_list, 'cluster': labels_list})
    return dataframe


def classify_clusters(out_dir, clusters_dir, n_threads):
    """
    Classify clusters with BAT

    :str out_dir: path to output files
    :str clusters_dir: path to directory of clustered fasta files (extension .fna)
    :int n_threads: number of threads
    :return: path to bin2classification files with names added
    """
    logging.info("Classifying clusters ... thank you for your patience")
    cat_out_path = os.path.join(out_dir, "cat-out")
    subprocess.check_output("{} bins -b {} -d {} -t {} --path_to_diamond {} -o {} -f 0.1 -s .fna"
                            " -n {}"
                            .format(CAT_EXEC, clusters_dir, CAT_DB, CAT_TAX,
                                    DIAMOND_EXEC, cat_out_path, n_threads), shell=True)
    orf2lca_file = "{}.ORF2LCA.txt".format(cat_out_path)
    orf2lca_names = "{}.ORF2LCA-names.txt".format(cat_out_path)
    bin2class_file = "{}.bin2classification.txt".format(cat_out_path)
    bin2class_names = "{}.bin2classification-names.txt".format(cat_out_path)
    subprocess.check_output("{} add_names -i {} -o {} -t {} --only_official"
                            .format(CAT_EXEC, orf2lca_file, orf2lca_names, CAT_TAX),
                            shell=True)
    subprocess.check_output("{} add_names -i {} -o {} -t {} --only_official"
                            .format(CAT_EXEC, bin2class_file, bin2class_names, CAT_TAX),
                            shell=True)
    return bin2class_names


def get_viral_bins(table_viral, df_clusters, bin2class_names_path):
    """
    Output information on bins containing viral sequences to infer phage interactions

    :str table_viral: path to .csv file with viral edges
    :pd.DataFrame df_clusters: Pandas df of node-cluster mapping
    :str bin2class_names_path: path to CAT bin2classification files with names added
    :return: Pandas df with viral nodes and their cluster information
    """
    logging.info("Placing viral sequences")
    df_viral = pd.read_csv(table_viral, dtype={'node': 'int32'}, sep='\t', index_col='node')
    df_clusters = df_clusters.astype({'node': 'int32'}).set_index('node')
    df_viral = df_viral.merge(df_clusters, on='node', how='left')
    if bin2class_names_path:
        df_bin2class = pd.read_csv(bin2class_names_path, sep='\t')
        df_bin2class['cluster'] = df_bin2class['# bin'].apply(lambda x: os.path.splitext(x)[0])
        df_bin2class = df_bin2class.drop(columns=['# bin', 'classification', 'reason'])
        df_bin2class = df_bin2class.astype({'cluster': 'int32'}).set_index('cluster')
        df_viral = df_viral.astype({'cluster': 'int32'}).set_index('cluster')
        df_viral = df_viral.merge(df_bin2class, on='cluster', how='left')
    return df_viral


def output_crispr_spacer_interactions(table_viral, df_clusters):
    """

    :param table_viral:
    :param df_clusters:
    :return:
    """
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    parser.add_argument(
        "input_graph", type=str,
        help="path to .gfa assembly graph")
    parser.add_argument(
        "--cat_pack", type=str, default=os.environ.get("CAT_PACK"),
        help="path to CAT_pack install")
    parser.add_argument(
        "--cat_db", type=str, default=os.environ.get("CAT_DB_DIR"),
        help="path to CAT database")
    parser.add_argument(
        "--cat_db_tax", type=str, default=os.environ.get("CAT_DB_TAX_DIR"),
        help="path to CAT taxonomy database")
    parser.add_argument(
        "--diamond_exec", type=str, default=os.environ.get("DIAMOND_EXEC"),
        help="path to DIAMOND executable")
    parser.add_argument(
        '--output-dir', type=str, default="./RHEA_results",
        help='output directory name [./RHEA_results]')
    parser.add_argument(
        '--avg-genome-len', type=int, default=4000000,
        help='average length of genomes in sample [4000000]')
    parser.add_argument(
        '--keep-cluster-fa', action="store_true",
        help='keep temp fasta file for each cluster')
    parser.add_argument(
        '--skip-classify', action="store_true",
        help='Skip CAT classification step')
    parser.add_argument(
        '--skip-clustering', action="store_true",
        help='Skip clustering')
    parser.add_argument(
        '--viral-table', type=str, default="",
        help='table of viral edges in graph with classifications')
    parser.add_argument(
        '--CRISPR', action="store_true",
        help='perform CRISPR spacer analysis')
    parser.add_argument(
        '--threads', type=int, default=3,
        help='threads utilized by DIAMOND [3]')
    args = parser.parse_args()

    # check CAT database and taxonomy are specified
    if not args.skip_classify:
        CAT_EXEC = os.path.join(args.cat_pack, "CAT")
        CAT_DB = args.cat_db
        CAT_TAX = args.cat_db_tax
        DIAMOND_EXEC = args.diamond_exec
        if not args.cat_pack:
            raise ValueError("CAT executable not specified. "
                             "Either 'export CAT_PACK=<path_to_CAT_pack>' or "
                             "utilize '--cat_pack' parameter.")
        if not CAT_DB:
            raise ValueError("CAT database not specified. "
                             "Either 'export CAT_DB_DIR=<path_to_database>' or "
                             "utilize '--cat_db' parameter.")
        if not CAT_TAX:
            raise ValueError("CAT database taxonomy not specified. "
                             "Either 'export CAT_DB_TAX_DIR=<path_to_database>' or "
                             "utilize '--cat_db_tax' parameter.")
        if not DIAMOND_EXEC:
            raise ValueError("DIAMOND executable not specified. "
                             "Either 'export DIAMOND_EXEC=<path_to_DIAMOND_exec>' or "
                             "utilize '--diamond_exec' parameter.")

    # create output directories
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_clusters_dir = os.path.join(args.output_dir, "clusters")
    if not os.path.exists(output_clusters_dir):
        os.makedirs(output_clusters_dir)

    # cluster graph
    if args.skip_clustering:
        graph_clusters_df = pd.read_csv(os.path.join(args.output_dir, "clusters.tsv"), sep='\t')
    else:
        weighted_graph, seqs = weight_graph(args.input_graph)
        graph_clusters_df = cluster_graph(weighted_graph, seqs, args.avg_genome_len,
                                          output_clusters_dir)
        graph_clusters_df.to_csv(os.path.join(args.output_dir, "clusters.tsv"), sep='\t', index=False)

    # classify clusters
    CLUSTER_CLASS = None
    if not args.skip_classify:
        CLUSTER_CLASS = classify_clusters(args.output_dir, output_clusters_dir,
                                          args.threads)
    # viral analysis
    if args.viral_table:
        viral_df = get_viral_bins(args.viral_table, graph_clusters_df, CLUSTER_CLASS)
        viral_df.to_csv(os.path.join(args.output_dir, "bins-with-viral.tsv"), sep='\t')
        if args.CRISPR:
            output_crispr_spacer_interactions(args.viral_table, graph_clusters_df)
