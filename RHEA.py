#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""

import os
import math
import argparse
import subprocess
import pandas as pd
import networkx as nx
from collections import defaultdict

import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import SpectralClustering


CAT_EXEC = "./exec/CAT"
DIAMOND_EXEC = "./exec/diamond"
__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'April 2023'

def weight_graph(graph_path):
    '''
    Parses .gfa file from MetaFlye output into a networkx graph. Nodes are weighted by length of sequence.
        edges are weighted by coverage.

    :param graph_path: (str) path to .gfa MetaFlye output
    :return: Weighted Networkx graph (G), dictionary of edge id to SeqRecord sequence (sequences)
    '''
    logging.info("Reading in graph: {}".format(graph_path))
    G = nx.Graph()
    sequences = {}
    with open(graph_path, "r") as f:

        for line in f:
            if line.startswith("L"):
                split_line = line.rstrip("\n").split("\t")
                node_A = split_line[1].split("_")[-1]
                node_B = split_line[3].split("_")[-1]
                RC = [cov.split(":")[-1] for cov in split_line if cov.startswith("RC")]
                weight = 0
                if RC:
                    weight = int(RC[0])
                # add weight to existing edge weight
                if not G.has_edge(node_A, node_B):
                    G.add_edge(node_A, node_B, weight=weight)
                else:
                    weight = G[node_A][node_B]["weight"] + weight
                    G.add_edge(node_A, node_B, weight=weight)

            if line.startswith("S"):
                split_line = line[:-1].split("\t")
                new_node = split_line[1].split("_")[-1]
                node_length = len(split_line[2])
                G.add_node(new_node, weight=node_length)
                sequences[new_node] = SeqRecord(Seq(split_line[2]), id=new_node, description="edge_{}".format(new_node))
    return G, sequences



def cluster_graph(G, sequences, genome_avg_len, output_path):
    '''
    Performs spectral clustering on weighted graph, where number of clusters decided such that each cluster is the
        length of one genome.

    :param G: weighted assembly graph. Nodes weighted by length of sequence. Edges weighted by coverage between nodes.
    :param sequences: dictionary of edge id to SeqRecord sequence
    :param genome_avg_len: expected average genome length for genomes in metagenome sample
    :param output_path: path to output clustered fasta files
    :return: Pandas df mapping graph nodes to assigned clusters
    '''
    logging.info("Clustering graphs, depending on the complexity of your graph, this may take several minutes;"
                 " output to: {}".format(output_path))
    nodes_list, labels_list = [], []
    total_bins = 0
    for cc in nx.connected_components(G):
        A = nx.to_numpy_array(G, nodelist=list(cc))
        length_sum = sum([G.nodes("weight")[node] for node in cc])
        n_clusters = math.ceil(length_sum/genome_avg_len)
        if n_clusters == 1:
            next_labels = [total_bins] * len(cc)
            cluster_seqs = [sequences.pop(node) for node in cc]
            SeqIO.write(cluster_seqs, os.path.join(output_path, "{}.fna".format(total_bins)), "fasta")
        else:
            sc = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
            sc.fit(A)
            next_labels = [i + total_bins for i in list(sc.labels_)]
            clusters_dict = defaultdict(list)
            for label, component_nodes in zip(next_labels, list(cc)):
                clusters_dict[label].append(component_nodes)
            for k, cluster_nodes in clusters_dict.items():
                cluster_seqs = [sequences.pop(node) for node in cluster_nodes]
                SeqIO.write(cluster_seqs, os.path.join(output_path, "{}.fa".format(k)), "fasta")
        nodes_list = nodes_list + list(cc)
        labels_list = labels_list + next_labels
        total_bins = total_bins + n_clusters
    df = pd.DataFrame({'node': nodes_list, 'cluster': labels_list})
    return df

def classify_clusters(out_dir, clusters_dir, cat_db, cat_db_tax, threads):
    '''
    Classify clusters with BAT

    :param out_dir: path to output files
    :param clusters_dir: path to directory of clustered fasta files (extension .fna)
    :param cat_db: path to CAT database
    :param cat_db_tax: path to CAT database taxonomy
    :return: path to bin2classifcation files with names added
    '''
    logging.info("Classifying clusters ... thank you for your patience")
    cat_out_path = os.path.join(out_dir, "cat-out")
    subprocess.check_output("{} bins -b {} -d {} -t {} --path_to_diamond {} -o {} -f 0.1 -s .fna -n {}"
                            .format(CAT_EXEC, clusters_dir, cat_db, cat_db_tax,
                                    DIAMOND_EXEC, cat_out_path, threads), shell=True)
    ORF2LCA_file = "{}.ORF2LCA.txt".format(cat_out_path)
    ORF2LCA_names = "{}.ORF2LCA-names.txt".format(cat_out_path)
    bin2class_file = "{}.bin2classification.txt".format(cat_out_path)
    bin2class_names = "{}.bin2classification-names.txt".format(cat_out_path)
    subprocess.check_output("{} add_names -i {} -o {} -t {} --only_official"
                            .format(CAT_EXEC, ORF2LCA_file, ORF2LCA_names, cat_db_tax), shell=True)
    subprocess.check_output("{} add_names -i {} -o {} -t {} --only_official"
                            .format(CAT_EXEC, bin2class_file, bin2class_names, cat_db_tax), shell=True)
    return bin2class_names


def get_viral_bins(table_viral, df_clusters, bin2class_names_path):
    '''
    Output information on bins containing viral sequences to infer phage interactions

    :param table_viral: path to .csv file with viral edges
    :param df_clusters: Pandas df of node-cluster mapping
    :param bin2class_names_path: path to CAT bin2classifcation files with names added
    :return: Pandas df with viral nodes and their cluster information
    '''
    logging.info("Detecting viral sequences")
    df_viral = pd.read_csv(table_viral, dtype={'node' : 'int32'})
    df_clusters = df_clusters.astype({'node' : 'int32'}).set_index('node')
    df_viral = df_viral.merge(df_clusters, on='node', how='left')
    if bin2class_names_path:
        df_bin2class = pd.read_csv(bin2class_names_path, sep='\t')
        df_bin2class['cluster'] = df_bin2class['# bin'].apply(lambda x: os.path.splitext(x)[0])
        df_bin2class = df_bin2class.drop(columns=['# bin', 'classification', 'reason'])
        df_bin2class = df_bin2class.astype({'cluster' : 'int32'}).set_index('cluster')
        df_viral = df_viral.astype({'cluster' : 'int32'}).set_index('cluster')
        df_viral = df_viral.merge(df_bin2class, on='cluster', how='left')
    return df_viral

def output_CRISPR_spacer_interactions(table_viral, df_clusters):
    return None



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    parser.add_argument(
        "input_graph", type=str,
        help="path to .gfa assembly graph")
    parser.add_argument(
        "--cat_db", type=str, default=os.environ.get("CAT_DB_DIR"),
        help="path to CAT database")
    parser.add_argument(
        "--cat_db_tax", type=str, default=os.environ.get("CAT_DB_TAX_DIR"),
        help="path to CAT taxonomy database")
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
        if not args.cat_db:
            raise ValueError("CAT database not specified. "
                         "Either 'export CAT_DB_DIR=<path_to_database>' or "
                         "utilize '--cat_db' parameter.")
        if not args.cat_db_tax:
            raise ValueError("CAT database taxonomy not specified. "
                         "Either 'export CAT_DB_TAX_DIR=<path_to_database>' or "
                         "utilize '--cat_db_tax' parameter.")

    # create output directories
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_clusters_dir = os.path.join(args.output_dir, "clusters")
    if not os.path.exists(output_clusters_dir):
        os.makedirs(output_clusters_dir)

    # cluster & classify graph
    weighted_graph, seqs = weight_graph(args.input_graph)
    graph_clusters_df = cluster_graph(weighted_graph, seqs, args.avg_genome_len, output_clusters_dir)
    graph_clusters_df.to_csv(os.path.join(args.output_dir, "clusters.tsv"), sep='\t', index=False)
    out_bin2class_names = None
    if not args.skip_classify:
        out_bin2class_names = classify_clusters(args.output_dir, output_clusters_dir, args.cat_db, args.cat_db_tax, args.threads)

    # viral analysis
    if args.viral_table:
        viral_df = get_viral_bins(args.viral_table, graph_clusters_df, out_bin2class_names)

        viral_df.to_csv(os.path.join(args.output_dir, "bins-with-viral.tsv"), sep='\t')
        if args.CRISPR:
            output_CRISPR_spacer_interactions(args.viral_table, graph_clusters_df)




