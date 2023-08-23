#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""

import os
import math
import shutil
import logging
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import SpectralClustering

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'May 2023'

TAX_RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
FASTA_EXT = "fna"


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
    nodes_data = []
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
                node_depth = int([cov.split(":")[-1] for cov in split_line if cov.startswith("dp")][0])
                new_node = split_line[1].split("_")[-1]
                node_length = len(split_line[2])
                graph.add_node(new_node, size=node_length, weight=node_depth)
                sequences[new_node] = SeqRecord(Seq(split_line[2]),
                                                id=new_node,
                                                description="edge_{}".format(new_node))
                nodes_data.append([new_node, node_length, node_depth])
    nodes_df = pd.DataFrame(nodes_data, columns=['node', 'node_length', 'node_depth'])
    return graph, sequences, nodes_df

def get_cc_stats(graph_cc, sum_length, n_nodes):

    """
    :param graph_cc: a graph that is a single connected component
    sum_length:
    n_nodes:

    :return: desired graph metrics
    """
    n_edges = len(graph_cc.edges)
    density = nx.density(graph_cc)
    cc_degree_dict = dict(graph_cc.degree())
    avg_degree = sum(cc_degree_dict.values()) / len(cc_degree_dict)
    highest_degree = max(cc_degree_dict.values())
    n_stars_d3 = sum(1 for value in cc_degree_dict.values() if value > 3)
    n_stars_d4 = sum(1 for value in cc_degree_dict.values() if value > 4)
    n_stars_d5 = sum(1 for value in cc_degree_dict.values() if value > 5)
    n_stars_d6 = sum(1 for value in cc_degree_dict.values() if value > 6)
    n_stars_d7 = sum(1 for value in cc_degree_dict.values() if value > 7)
    n_stars_d10 = sum(1 for value in cc_degree_dict.values() if value > 10)

    clustering_coef = nx.average_clustering(graph_cc)
    n_self_loop = nx.number_of_selfloops(graph_cc)
    nodes_bps = nx.get_node_attributes(graph_cc, 'size').values()
    avg_node_bp =  sum(nodes_bps) / len(nodes_bps)
    median_node_bp = np.median(list(nodes_bps))
    edge_weights = nx.get_edge_attributes(graph_cc, 'weight').values()
    avg_edge_weight, median_edge_weight = np.nan, np.nan
    if len(edge_weights) > 0 :
        avg_edge_weight = sum(edge_weights) / len(edge_weights)
        median_edge_weight = np.median(list(edge_weights))
    betweenness_dict_cc = nx.betweenness_centrality(graph_cc)
    avg_betweenness = sum(betweenness_dict_cc.values()) / len(betweenness_dict_cc)
    # nx.approximate_current_flow_betweenness_centrality(graph_cc)
    triangles_dict_cc = nx.triangles(graph_cc)
    triangles_count = sum(triangles_dict_cc.values()) / 3
    return {"n nodes": n_nodes, "n edges": n_edges, "total length (bp)": sum_length,
            "mean node length": avg_node_bp, "median node length": median_node_bp,
            "mean edge weight": avg_edge_weight, "median edge weight": median_edge_weight,
            "density": density, "mean degree": avg_degree, "highest degree": highest_degree,
            "n self loops": n_self_loop, "clustering coefficient": clustering_coef,
            "mean betweenness": avg_betweenness, "n triangles": triangles_count,
            "n stars 3": n_stars_d3, "n stars 4": n_stars_d4, "n stars 5": n_stars_d5,
            "n stars 6": n_stars_d6, "n stars 7": n_stars_d7, "n stars 10": n_stars_d10}, cc_degree_dict, betweenness_dict_cc, triangles_dict_cc


def cluster_graph(graph, sequences, genome_avg_len, genome_min_len,
                  nodes_df, output_path_clusters, output_path_fragments,
                  skip_output_fasta):
    """
    Performs spectral clustering on weighted graph, where number of clusters decided such that
    each cluster is the length of one genome.

    :networkx.G graph: weighted assembly graph. Nodes size length of sequence and weighted by coverage.
            Edges weighted by coverage between nodes.
    :dict{str:SeqRecord} sequences: dictionary of edge id to SeqRecord sequence
    :int genome_avg_len: expected average genome length for genomes in metagenome sample
    :int genome_min_len: minimum genome length for genomes in metagenome sample. Used as a
            threshold for clusters vs fragments
    :list of lists nodes_data:
    :str output_path_clusters: path to output clustered fasta files
    :str output_path_fragments: path to output fragments fasta files
    :return: Pandas df mapping graph nodes to assigned clusters
    """
    logging.info("Clustering graphs, depending on the complexity of your graph, "
                 "this may take several minutes; output to: %s", output_path_clusters)
    all_nodes, cluster_labels = [], []
    cut_edges_nodes, cut_edges_weights, cut_edges_clusters = [], [], []
    degree_dict, betweenness_dict, triangles_dict = {}, {}, {}
    fragment_nodes, fragment_labels = [], []
    cluster_stats = {}  # track number of nodes and cumulative node weight for each cluster
    total_bins = 0
    for connected_component in nx.connected_components(graph):
        adj_matrix = nx.to_numpy_array(graph, nodelist=list(connected_component))
        length_sum = sum([int(graph.nodes("size")[node]) for node in connected_component])
        n_clusters = math.ceil(length_sum / genome_avg_len)
        if n_clusters == 1 or len(connected_component) < 3:  # connected component is one cluster or fragment
            next_labels = [total_bins] * len(connected_component)
            if sequences:
                cluster_seqs = [sequences.pop(node) for node in connected_component]
            if length_sum < genome_min_len:  # add to fragments list if len_sum under threshold
                if sequences and not skip_output_fasta:
                    SeqIO.write(cluster_seqs, os.path.join(output_path_fragments, "{}.{}".format(total_bins,FASTA_EXT)),
                                "fasta")
                fragment_nodes = fragment_nodes + list(connected_component)
                fragment_labels = fragment_labels + next_labels
            else:  # else add to clusters
                if sequences and not skip_output_fasta:
                    SeqIO.write(cluster_seqs, os.path.join(output_path_clusters, "{}.{}".format(total_bins,FASTA_EXT)),
                                "fasta")
                all_nodes = all_nodes + list(connected_component)
                cluster_labels = cluster_labels + next_labels
                cuts = [[] for _ in range(len(next_labels))]
                cut_edges_nodes = cut_edges_nodes + cuts
                cut_edges_clusters = cut_edges_clusters + cuts
                cut_edges_weights = cut_edges_weights + cuts
                cluster_stats[total_bins], cc_degree_dict, cc_betweenness_dict, cc_triangles_dict = get_cc_stats(graph.subgraph(connected_component),
                                                     length_sum, len(connected_component))
                degree_dict.update(cc_degree_dict)
                betweenness_dict.update(cc_betweenness_dict)
                triangles_dict.update(cc_triangles_dict)



        else:  # if cc has more than one cluster
            spec_clust = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
            spec_clust.fit(adj_matrix)
            next_labels = [i + total_bins for i in list(spec_clust.labels_)]
            clusters_dict = defaultdict(list)
            for label, component_nodes in zip(next_labels, list(connected_component)):
                clusters_dict[label].append(component_nodes)
            for k, cluster_nodes in clusters_dict.items():
                length_cluster = sum([graph.nodes("size")[node] for node in cluster_nodes])
                cluster_stats[k], cc_degree_dict, cc_betweenness_dict, cc_triangles_dict = get_cc_stats(graph.subgraph(cluster_nodes),
                                                         length_cluster, len(cluster_nodes))
                degree_dict.update(cc_degree_dict)
                betweenness_dict.update(cc_betweenness_dict)
                triangles_dict.update(cc_triangles_dict)
                if sequences and not skip_output_fasta:
                    cluster_seqs = [sequences.pop(node) for node in cluster_nodes]
                    SeqIO.write(cluster_seqs, os.path.join(output_path_clusters, "{}.{}".format(k, FASTA_EXT)), "fasta")
            all_nodes = all_nodes + list(connected_component)
            cluster_labels = cluster_labels + next_labels

            ## adj_matrix after spectral clustering
            n_nodes = len(list(connected_component))
            clustered_adj_matrix = np.zeros((n_nodes, n_nodes), dtype=int)
            for i in range(n_nodes):
                for j in range(n_nodes):
                    if cluster_labels[i] == cluster_labels[j]:
                        clustered_adj_matrix[i, j] = adj_matrix[i, j]
            cut_nodes_matrix = adj_matrix - clustered_adj_matrix

            ## get information on cut edges
            for i in range(n_nodes):
                cut_indexes = [int(val) for val in np.where(cut_nodes_matrix[i])[0]]
                cut_edges_nodes.append([int(val) for val in [list(connected_component)[cut_index] for cut_index in cut_indexes]])
                cut_edges_clusters.append([next_labels[cut_index] for cut_index in cut_indexes])
                cut_edges_weights.append([int(val) for val in cut_nodes_matrix[i][cut_nodes_matrix[i] != 0]])
        total_bins = total_bins + n_clusters

    clusters_dataframe = pd.DataFrame({"node": all_nodes,
                                       "cluster": cluster_labels})
    clusters_dataframe = clusters_dataframe.merge(nodes_df, on="node", how="left")
    clusters_dataframe["cut edges"] = cut_edges_nodes
    clusters_dataframe["cut edges weights"] = cut_edges_weights
    clusters_dataframe["cut edges clusters"] = cut_edges_clusters
    clusters_dataframe = clusters_dataframe.sort_values(['cluster', 'node'])
    clusters_dataframe.set_index("node", inplace=True)
    clusters_dataframe['degree'] = clusters_dataframe.index.map(degree_dict)
    clusters_dataframe['betweenness'] = clusters_dataframe.index.map(betweenness_dict)
    clusters_dataframe['triangles'] = clusters_dataframe.index.map(triangles_dict)

    fragments_dataframe = pd.DataFrame({"node": fragment_nodes,
                                        "cluster": fragment_labels})
    fragments_dataframe = fragments_dataframe.merge(nodes_df, on="node", how="left")
    cluster_info_df = pd.DataFrame(cluster_stats.values(), index=cluster_stats.keys())
    return clusters_dataframe, \
           fragments_dataframe.sort_values(['cluster', 'node']), \
           cluster_info_df.reset_index().rename(columns={'index': 'cluster'}).sort_values('cluster')


def classify_clusters(out_dir, clusters_dir, group_type, n_threads):
    """
    Classify clusters with BAT

    :str out_dir: path to output files
    :str clusters_dir: path to directory of clustered fasta files (extension .fna)
    :str group_type: denote if clusters or fragments
    :int n_threads: number of threads
    :return: path to bin2classification files with names added
    """
    logging.info("Classifying %s ... thank you for your patience", group_type)
    cat_out_path = os.path.join(out_dir, "cat-out-{}".format(group_type))
    subprocess.check_output("{} bins -b {} -d {} -t {} --path_to_diamond {} -o {} -f 0.1 -n {}"
                            " -s {}"
                            .format(CAT_EXEC, clusters_dir, CAT_DB, CAT_TAX,
                                    DIAMOND_EXEC, cat_out_path, n_threads, FASTA_EXT), shell=True)
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


def merge_classification_into_clusters_info(df_cluster_info, clusters_class_path):
    cat_df = pd.read_csv(clusters_class_path, sep='\t')
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    merge_on = "cluster"
    cat_df = cat_df.fillna("NaN")
    cat_df["count"] = 1
    cat_df["cluster"] = cat_df["# bin"].apply(lambda x: x.split(".")[0]).drop(columns="# bin")
    for col in ranks:
        cat_df[col] = cat_df[col].apply(lambda x: x.split(":")[0])
    agg_dict = {col: lambda x: ', '.join(x) if x.nunique() > 1 else x.iloc[0] for col in cat_df.columns if col != merge_on}
    agg_dict['count'] = 'sum'  # Add a new aggregation function to count rows
    cat_df_merged = cat_df.groupby(merge_on).agg(agg_dict).reset_index()
    for col in ranks:
        df_cluster_info[col] = cat_df_merged[col]
    df_cluster_info['count'] = cat_df_merged['count']
    return df_cluster_info


def get_viral_bins(table_viral, df_clusters, bin2class_names_path, group_type):
    """
    Output information on bins containing viral sequences to infer phage interactions

    :str table_viral: path to .csv file with viral edges
    :pd.DataFrame df_clusters: Pandas df of node-cluster mapping
    :str bin2class_names_path: path to CAT bin2classification files with names added
    :str group_type: denote if clusters or fragments
    :return: Pandas df with viral nodes and their cluster information
    """

    logging.info("Placing viral sequences in {}".format(group_type))
    df_viral = pd.read_csv(table_viral, dtype={"node": "int32"}, sep="\t", index_col="node")
    df_clusters = df_clusters.astype({"node": "int32"}).set_index("node")
    df_viral = df_viral.merge(df_clusters, on="node", how="inner").reset_index()
    if bin2class_names_path:
        df_bin2class = pd.read_csv(bin2class_names_path, sep="\t")
        df_bin2class["cluster"] = df_bin2class["# bin"].apply(lambda x: os.path.splitext(x)[0])
        df_bin2class = df_bin2class.astype({"cluster": "int32"})
        df_bin2class_flat = df_bin2class[["cluster"] + TAX_RANKS]
        df_bin2class_flat = df_bin2class_flat.groupby("cluster", as_index=False)[TAX_RANKS].agg(lambda x: list(x))
        df_viral = df_viral.astype({"cluster": "int32"})
        df_viral = df_viral.merge(df_bin2class_flat, on="cluster", how="left")
    return df_viral


def output_crispr(dir_clusters, dir_minced_output, minced_exec, threads):
    """ Run MinCED to detect all CRISPR in RHEA clusters

    :param dir_clusters
    :param dir_minced_output
    :param minced_path
    :return:
    """
    # run MinCED on each cluster
    logging.info("Running MinCED on each cluster; output to: %s", dir_minced_output)
    for file in os.listdir(dir_clusters):
        if file.endswith(FASTA_EXT):
            file_stem = Path(file).stem
            output_file = os.path.join(dir_minced_output, "{}.minced".format(file_stem))
            subprocess.check_output("{} {} {}".format(minced_exec, os.path.join(dir_clusters, file), output_file), shell=True)
    return None


def output_crispr_spacer_interactions(dir_fragments, dir_minced_output, dir_spacepharer_output,
                                      spacepharer_exec, threads):
    """ Run spacepharer between the CRISPR spacers detected in clusters and the RHEA fragments

    :return:
    """
    # run spacepharer between CRISPR in clusters and fragments
    logging.info("Running spacepharer; output to: %s", dir_spacepharer_output)
    fragments_all = os.path.join(dir_fragments, "*.{}".format(FASTA_EXT))
    minced_all = os.path.join(dir_minced_output, "*.{}".format("minced"))
    curr_dir = os.getcwd()
    os.chdir(dir_spacepharer_output)

    subprocess.check_output("{} createsetdb {} fragmentTargetSetDB spacepharer_tmp --threads {}"
                            .format(spacepharer_exec, fragments_all, threads), shell=True)
    subprocess.check_output("{} createsetdb {} fragmentControlSetDB spacepharer_tmp --reverse-fragments 1 --threads {}"
                            .format(spacepharer_exec, fragments_all, threads), shell=True)
    subprocess.check_output("{} parsespacer {} queryDB --threads {}"
                            .format(spacepharer_exec, minced_all, threads), shell=True)
    subprocess.check_output("{} createsetdb queryDB querySetDB spacepharer_tmp --extractorf-spacer 1 --threads {}"
                            .format(spacepharer_exec, threads), shell=True)
    subprocess.check_output("{} predictmatch querySetDB fragmentTargetSetDB fragmentControlSetDB "
                            "spacepharer_results.tsv spacepharer_tmp --threads {}"
                            .format(spacepharer_exec, threads), shell=True)
    results_path = os.path.join(os.getcwd(), "spacepharer_results.tsv")
    os.chdir(curr_dir)
    return results_path



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
        "--minced_exec", type=str, default="minced",
        help="path to MinCED exectuable")
    parser.add_argument(
        "--spacepharer_exec", type=str, default="spacepharer",
        help="path to spacepharer exectuable")
    parser.add_argument(
        '--output-dir', type=str, default=os.path.join(os.getcwd(), "RHEA_results"),
        help='output directory name [./RHEA_results]')
    parser.add_argument(
        '--avg-genome-len', type=int, default=4000000,
        help='average length of genomes in sample [4000000]')
    parser.add_argument(
        '--min-genome-len', type=int, default=400000,
        help='min length of genomes in sample [400000]')
    parser.add_argument(
        '--keep-cluster-fa', action="store_true",
        help='keep temp fasta file for each cluster')
    parser.add_argument(
        '--input-classifications', nargs=2, metavar=('clusters_class', 'fragments_class'),
        help='input classifications if skipping CAT classification step')
    parser.add_argument(
        '--skip-classify', action="store_true",
        help='Skip CAT classification step')
    parser.add_argument(
        '--skip-clustering', action="store_true",
        help='Skip clustering')
    parser.add_argument(
        '--skip-output-fasta', action="store_true",
        help='Skip additional fasta output for each cluster')
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

    # if CRISPR detection is on, check skip output fasta is off and executables
    if args.CRISPR:
        if args.skip_output_fasta:
            raise ValueError("skip-output-fasta must be set to false "
                         "in CRISPR detection mode.")
        # check minced and spacepharer exectuables are active


    # create output directories
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    if not os.path.isabs(args.output_dir):
        args.output_dir = os.path.normpath(os.path.join(os.getcwd(), args.output_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    output_clusters_dir = os.path.join(args.output_dir, "clusters")
    output_fragments_dir = os.path.join(output_clusters_dir, "fragments")
    output_minced_dir = os.path.join(args.output_dir, "minced")
    output_spacepharer_dir = os.path.join(args.output_dir, "spacepharer")
    if os.path.exists(output_spacepharer_dir):
        shutil.rmtree(output_spacepharer_dir)
        os.makedirs(output_spacepharer_dir)


    # cluster graph
    if args.skip_clustering:
        # read in existing clusters mapping file
        graph_clusters_df = pd.read_csv(os.path.join(args.output_dir, "clusters.tsv"), sep='\t')
        graph_fragments_df = pd.read_csv(os.path.join(args.output_dir, "fragments.tsv"), sep='\t')
    else:
        # create directories for clusters and fragments
        if not os.path.exists(output_clusters_dir):
            os.makedirs(output_clusters_dir)
        if not os.path.exists(output_fragments_dir):
            os.makedirs(output_fragments_dir)
        # create and cluster RHEA graph
        if os.path.splitext(args.input_graph)[1] == '.graphml':
            seqs = None
            df_nodes = pd.DataFrame([], columns=['node'])
            weighted_graph = nx.read_graphml(args.input_graph)

        else:
            weighted_graph, seqs, df_nodes = weight_graph(args.input_graph)
        graph_clusters_df, graph_fragments_df, cluster_info_df = cluster_graph(weighted_graph, seqs,
                                                              args.avg_genome_len,
                                                              args.min_genome_len,
                                                              df_nodes,
                                                              output_clusters_dir,
                                                              output_fragments_dir,
                                                              args.skip_output_fasta)
        graph_clusters_df.to_csv(os.path.join(args.output_dir, "clusters.tsv"), sep='\t', index=False)
        graph_fragments_df.to_csv(os.path.join(args.output_dir, "fragments.tsv"), sep='\t', index=False)
        cluster_info_df.to_csv(os.path.join(args.output_dir, "clusters_info.tsv"), sep='\t', index=False)


    # classify clusters
    if not args.skip_classify:
        CLUSTER_CLASSIFICATIONS = classify_clusters(args.output_dir, output_clusters_dir,
                                                    "clusters", args.threads)
        FRAGMENT_CLASSIFICATIONS = classify_clusters(args.output_dir, output_fragments_dir,
                                                     "fragments", args.threads)
    elif args.input_classifications:
        clusters_class, fragments_class = args.input_classifications
        CLUSTER_CLASSIFICATIONS = clusters_class
        FRAGMENT_CLASSIFICATIONS = fragments_class
    else:
        CLUSTER_CLASSIFICATIONS, FRAGMENT_CLASSIFICATIONS = None, None

    if CLUSTER_CLASSIFICATIONS and not args.skip_clustering:
        cluster_info_df = merge_classification_into_clusters_info(cluster_info_df, CLUSTER_CLASSIFICATIONS)
        cluster_info_df.to_csv(os.path.join(args.output_dir, "clusters_info.tsv"), sep='\t', index=False)



    # viral analysis
    if args.viral_table:
        viral_df_clusters = get_viral_bins(args.viral_table, graph_clusters_df,
                                           CLUSTER_CLASSIFICATIONS, "cluster")
        viral_df_clusters.to_csv(os.path.join(args.output_dir, "clusters-with-viral.tsv"), sep='\t', index=False)
        viral_df_fragments = get_viral_bins(args.viral_table, graph_fragments_df,
                                            FRAGMENT_CLASSIFICATIONS, "fragment")
        viral_df_fragments.to_csv(os.path.join(args.output_dir, "fragments-with-viral.tsv"), sep='\t', index=False)
    if args.CRISPR:
        if not os.path.exists(output_minced_dir):
            os.makedirs(output_minced_dir)
        if not os.path.exists(output_spacepharer_dir):
            os.makedirs(output_spacepharer_dir)
        output_crispr(output_clusters_dir, output_minced_dir, args.minced_exec, args.threads)
        spacepharer_results_path = output_crispr_spacer_interactions(output_fragments_dir, output_minced_dir,
                                          output_spacepharer_dir, args.spacepharer_exec, args.threads)
        spacepharer_destination_path = os.path.join(args.output_dir, os.path.basename(spacepharer_results_path))
        os.rename(spacepharer_results_path, spacepharer_destination_path)

