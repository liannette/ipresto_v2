import os
import random
import time
from collections import Counter, OrderedDict
from functools import partial
from itertools import combinations, islice
from multiprocessing import Pool

import networkx as nx
from sympy import binomial as ncr

from ipresto.preprocess.utils import (
    format_cluster_to_string,
    read_clusters,
    write_gene_counts,
    write_representatives,
)


def generate_adjacent_domain_pairs(domains):
    """
    Generate a set of unique pairs of adjacent domains, excluding pairs where one domain is '-'.

    Args:
        domains (list of str): A list of domain names, including '-' to represent an empty gene.

    Returns:
        set: A set of tuples, each containing a pair of adjacent domain names in sorted order.
             Pairs that include a domain separated by '-' are excluded.
    """
    return {
        tuple(sorted(pair))
        for pair in zip(domains[:-1], domains[1:])
        if "-" not in pair
    }


def calc_adj_index(domains1, domains2):
    """
    Calculate the adjacency index between two clusters.

    The adjacency index is a measure of similarity between two clusters based on
    the adjacency of their domains. Two domains are considered adjacent if they
    are consecutive in the cluster and are not separated by a '-' character.

    Parameters:
    doms1 (list of str): The first cluster, represented as a list of domain strings.
    doms2 (list of str): The second cluster, represented as a list of domain strings.

    Returns:
    float: The adjacency index, a value between 0.0 and 1.0, where 0.0 indicates no
           adjacency and 1.0 indicates perfect adjacency between the two clusters.

    Note:
    - If either cluster has no adjacent domain pairs (i.e., all domains are separated
      by '-'), the function returns 0.0.


    Returns the adjacency index between two clusters

    doms1, doms2: list of str, domainlist of a cluster

    If there is an empty gene between two domains these two domains are not
        adjacent
    """
    domain_pairs1 = generate_adjacent_domain_pairs(domains1)
    domain_pairs2 = generate_adjacent_domain_pairs(domains2)

    # If either cluster has no domain pairs, the adjacency index is 0.0
    if not domain_pairs1 or not domain_pairs2:
        return 0.0

    # Calculate adjacency index as ratio of intersection to union of domain pairs
    adj_index = len(domain_pairs1 & domain_pairs2) / len(domain_pairs1 | domain_pairs2)
    return adj_index


def is_contained(doms1, doms2):
    """
    Check if all domains from one cluster are contained in the other cluster.

    Parameters:
    doms1 (list of str): Domain list of the first cluster.
    doms2 (list of str): Domain list of the second cluster.

    Returns:
    bool: True if all domains from one cluster are contained in the other, False otherwise.
    """
    all_in_doms2 = all(domain in doms2 for domain in doms1 if domain != "-")
    all_in_doms1 = all(domain in doms1 for domain in doms2 if domain != "-")
    return all_in_doms2 or all_in_doms1


def extract_domains(tokenised_genes):
    """
    Extracts the domains from a list of tokenised genes.

    Args:
        tokenised_genes (list): A list of tokenised genes.

    Returns:
        list: A list of domains.
    """
    return [domain for gene in tokenised_genes for domain in gene]


def generate_edge(bgc_pair, adj_cutoff):
    """
    Calculate similarity scores between two bgcs and return if above cutoff

    pair: tuple of 2 bgcs
    cutoff: float
    A tuple is returned that can be read as an edge by nx.Graph.add_edges_from
    """
    bgc_name1, genes1 = bgc_pair[0]
    bgc_name2, genes2 = bgc_pair[1]

    domains1 = extract_domains(genes1)
    domains2 = extract_domains(genes2)

    contained = is_contained(domains1, domains2)
    adj_index = calc_adj_index(domains1, domains2)

    if contained or adj_index > adj_cutoff:
        return bgc_name1, bgc_name2, adj_index, contained


def generate_edges(dom_dict, cutoff, cores, temp_file, verbose):
    """Returns a pair of clusters in a tuple if ai/contained above cutoff

    dom_dict: dict {clus1:[domains]}, clusters linked to domains
    cutoff: float, between 0-1, when clusters are similar
    cores: int, amount of cores used for calculation

    returns a generator

    Note: temp file storing the edges so they are not in memory and passed in pool
    """
    if verbose:
        print("\nGenerating similarity scores")

    clusters = dom_dict.items()
    pairs = combinations(clusters, 2)
    slice_size = int(ncr(25000, 2))
    tot_size = ncr(len(clusters), 2)
    slce = islice(pairs, slice_size)
    chunk_num = int(tot_size / slice_size) + 1
    tloop = time.time()
    # update tempfile with increments of slice_size
    for i in range(chunk_num):
        if i == chunk_num - 1:
            # get chunksize of remainder
            chnksize = int(
                ((tot_size / slice_size % 1 * slice_size) / (cores * 20)) + 1
            )
            if chnksize < 5:
                chnksize = 5
        else:
            # the default used by map divided by 5
            chnksize = int((slice_size / (cores * 20)) + 1)
        pool = Pool(cores, maxtasksperchild=10)
        edges_slce = pool.imap(
            partial(generate_edge, adj_cutoff=cutoff), slce, chunksize=chnksize
        )
        pool.close()
        pool.join()
        # write to file
        with open(temp_file, "a") as tempf:
            for line in edges_slce:
                if line:
                    tempf.write("{}\n".format("\t".join(map(str, line))))
        slce = islice(pairs, slice_size)
        del (edges_slce, pool)
        if verbose and i == 0:
            t = (time.time() - tloop) * chunk_num
            # based on one loop
            t_str = "  it will take around {}h{}m{}s".format(
                int(t / 3600), int(t % 3600 / 60), int(t % 3600 % 60)
            )
            print(t_str)


def generate_graph(edges, verbose):
    """Returns a networkx graph

    edges: list/generator of tuples, (pair1,pair2,{attributes})
    """
    g = nx.Graph()
    g.add_edges_from(edges)
    if verbose:
        print("\nGenerated graph with:")
        print(" {} nodes".format(g.number_of_nodes()))
        print(" {} edges".format(g.number_of_edges()))
    return g


def read_edges(file_path):
    """Yields edges from temp file

    file_path: str
    """
    with open(file_path, "r") as inf:
        for line in inf:
            line = line.strip("\n").split("\t")
            cont = line[-1] == "True"
            tup = (line[0], line[1], {"ai": float(line[2]), "contained": cont})
            yield tup


def find_representatives(cliques, d_l_dict, graph):
    """
    Returns {representative:[clique]} based on cluster/bgc with most domains in clique

    cliques: list of lists of strings, cliques of clusters
    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    graph: networkx graph structure of the cliques
    The longest cluster is chosen (most domains). If there are multiple
        longest clusters then the cluster with the least connections is
        chosen (to preserve most information).
    """
    representative_bgcs = OrderedDict()
    redundant_bgcs = set()

    for cliq in cliques:
        # Filter out already processed BGCs
        cliq = [bgc for bgc in cliq if bgc not in redundant_bgcs]
        if not cliq:
            continue

        # Find the BGC with the maximum number of domains
        domlist = [(bgc, d_l_dict[bgc]) for bgc in cliq]
        maxdoml = max(doml for _, doml in domlist)
        clus_maxlen = [bgc for bgc, doml in domlist if doml == maxdoml]

        # If there are multiple, choose the one with the minimum degree
        if len(clus_maxlen) > 1:
            min_degr = min(graph.degree(bgc) for bgc in clus_maxlen)
            rep = random.choice(
                [bgc for bgc in clus_maxlen if graph.degree(bgc) == min_degr]
            )
        else:
            rep = clus_maxlen[0]

        # Update the representative BGCs
        if rep not in representative_bgcs:
            representative_bgcs[rep] = set()
        representative_bgcs[rep].update(cliq)

        # Mark the remaining BGCs as processed
        cliq.remove(rep)
        redundant_bgcs.update(cliq)

    return representative_bgcs


def find_all_representatives(d_l_dict, g):
    """Iterates find_representatives until there are no similar bgcs

    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    g: networkx graph structure containing the cliques
    all_reps_dict: dict of {representative:[represented]}
    """
    print("\nFiltering out similar bgcs.")

    all_reps_dict = {}
    subg = g.subgraph(g.nodes)
    i = 1
    while subg.number_of_edges() != 0:
        print(
            f"  iteration {i}, edges (similarities between bgcs) left: {subg.number_of_edges()}"
        )
        cliqs = nx.algorithms.clique.find_cliques(subg)
        # make reproducible by making the cliqs have the same order every time
        # sort first each cliq alphabetically, then cliqs alphabetically,
        # then on length, so longest are first and order is the same
        cliqs = sorted(sorted(cl) for cl in cliqs if len(cl) > 1)
        cliqs.sort(key=len, reverse=True)
        reps_dict = find_representatives(cliqs, d_l_dict, subg)
        subg = subg.subgraph(reps_dict.keys())
        # merge reps_dict with all_reps_dict
        for key, vals in reps_dict.items():
            if key not in all_reps_dict:
                all_reps_dict[key] = vals
            else:
                # merge represented clusters in a new representative
                newvals = []
                for old_rep in vals:
                    # if statement for bgcs already represented by this
                    # representative and thus no longer in all_reps_dict
                    if old_rep in all_reps_dict.keys():
                        newv = [v for v in all_reps_dict[old_rep]]
                        newvals += newv
                        del all_reps_dict[old_rep]
                all_reps_dict[key] = set(newvals)
        i += 1
    return all_reps_dict


def calculate_domain_lengths(tokenised_clusters):
    """
    Calculates the number of domains for each BGC, excluding empty genes ("-").

    Args:
        tokenised_clusters (dict): A dictionary where keys are BGCs and values are lists of domains.

    Returns:
        dict: A dictionary where keys are BGCs and values are the number of domains.
    """
    return {
        bgc: sum(len(g) for g in genes if g != ("-",))
        for bgc, genes in tokenised_clusters.items()
    }


def get_representatives(clusters, graph):
    """
    Identifies representative clusters from the given tokenised clusters.

    Args:
        clusters (dict): A dictionary where keys are BGCs and values are lists of domains.
        graph (nx.Graph): A graph representing the similarity between clusters.
        sim_cutoff (float): Similarity cutoff for generating edges.
        cores (int): Number of cores to use for parallel processing.

    Returns:
        dict: A dictionary where keys are representative BGCs and values are lists of associated BGCs.
    """
    domains_per_bgc = calculate_domain_lengths(clusters)
    representative_bgcs = find_all_representatives(domains_per_bgc, graph)
    unique_bgcs = {bgc: [bgc] for bgc in clusters if bgc not in graph}
    representative_bgcs.update(unique_bgcs)
    return representative_bgcs


def filter_similar_clusters(
    in_file_path: str,
    out_file_path: str,
    counts_file_path: str,
    representatives_file_path: str,
    edges_file_path: str,
    sim_cutoff: float,
    cores: int,
    verbose: bool,
) -> str:
    """Removes clusters based on similarity and writes the filtered clusters to a file.

    This function reads a file containing tokenised clusters, filters out clusters based on a
    similarity cutoff, and writes the filtered clusters to a new file. If redundancy
    filtering is disabled or if the filtered file already exists, the function skips
    the filtering process.
    """
    if verbose:
        print(f"\nPerforming similarity filtering on {in_file_path}.")

    # Read clusters from input file
    clusters = read_clusters(in_file_path)
    # Generate edges between clusters based on similarity
    generate_edges(clusters, sim_cutoff, cores, edges_file_path, verbose)
    # Generate a graph from the edges
    graph = generate_graph(read_edges(edges_file_path), verbose)

    # Find representative clusters
    representative_bgcs = get_representatives(clusters, graph)
    write_representatives(representative_bgcs, representatives_file_path)

    # Count the occurence of each gene
    gene_counter = Counter()
    for cluster_id in representative_bgcs:
        genes = clusters[cluster_id]
        gene_counter.update(genes)

    # Write the results to the output files
    write_gene_counts(gene_counter, counts_file_path)
    with open(out_file_path, "w") as outfile:
        for cluster_id in representative_bgcs:
            genes = clusters[cluster_id]
            outfile.write(format_cluster_to_string(cluster_id, genes))

    # Clean up temporary files
    os.remove(edges_file_path)

    if verbose:
        print("\nSimilarity filtering complete.")
        print(
            f"Selected {len(representative_bgcs)} representatives for {len(clusters)} clusters"
        )
        print(f"Representative clusters written to {out_file_path}")
