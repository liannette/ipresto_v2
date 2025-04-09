import re

from ipresto.preprocess.utils import (
    count_non_emtpy_genes,
    parse_cluster_line,
    format_cluster_to_string,
    write_gene_counts,
    read_txt,
)
from collections import Counter
from typing import Tuple, Set
from multiprocessing import Pool
from functools import partial


def extract_domain_base_name(domain: str) -> str:
    """Extract the base name of a domain.

    This function extracts the base name of a protein_domain,
    excluding any '_c' suffix of subPfams.

    Args:
        domain (str): The domain name.

    Returns:
        str: The base name of the domain.
    """
    return re.sub(r"_c\d+$", "", domain)


def filter_domains_in_gene(gene: Tuple[str], include_domains: Set[str]) -> Tuple[str]:
    """Filter domains in a gene based on the include_list.

    This function iterates through the domains in a given gene and filters
    them based on whether their base name (excluding any '_c' suffix) is
    present in the include_list. If a domain matches the criteria, it is
    included in the filtered_gene list.

    Args:
        gene (list of str): A list of domain names in the gene.
        include_domains (set of str): A set of domain base names to include.

    Returns:
        list of str: A list of filtered domain names. If no domains match
            the criteria, a list containing a single tuple with a dash ('-')
            is returned.
    """
    filtered_gene = []
    for domain in gene:
        domain_name = extract_domain_base_name(domain)
        if domain_name in include_domains:
            filtered_gene.append(domain)
    return tuple(filtered_gene) if filtered_gene else ("-",)


def process_cluster(cluster, include_domains, min_genes, verbose):
    cluster_id, genes = cluster
    filtered_genes = [filter_domains_in_gene(gene, include_domains) for gene in genes]
    n_genes_with_domains = count_non_emtpy_genes(filtered_genes)
    if n_genes_with_domains < min_genes:
        if verbose:
            print(f"  excluding {cluster_id}, only {n_genes_with_domains} genes with domain hits (min {min_genes})")
        return None
    else:
        return cluster_id, filtered_genes


def perform_domain_filtering(
    in_file_path: str,
    domain_filtering_file_path: str,
    out_file_path: str,
    counts_file_path: str,
    min_genes: int,
    cores: int,
    verbose: bool,
) -> str:
    """
    Wrapper for domain filtering of clusters.

    Args:
        in_file_path (str): Path to the input file containing tokenised clusters.
        domain_filtering_file_path (str): Path to the file containing the list of domains to include.
        out_file_path (str): Path to the output file for writing the domain-filtered clusters.
        counts_file_path (str): Path to the file for writing gene counts.
        min_genes (int): Minimum number of non-empty genes required in a cluster.
        verbose (bool): If True, print verbose output.
    """
    if verbose:
        print(f"\nPerforming domain filtering on {in_file_path}")
        print(f"Removing protein domains not listed in {domain_filtering_file_path}")

    # Read the input files
    include_domains = set(read_txt(domain_filtering_file_path))
    with open(in_file_path, "r") as infile:
        clusters = [parse_cluster_line(line) for line in infile.readlines()]

    # Process each cluster in parallel
    with Pool(cores, maxtasksperchild=1000) as pool:
        process_func = partial(
            process_cluster,
            include_domains=include_domains,
            min_genes=min_genes,
            verbose=verbose,
        )
        results = pool.map(process_func, clusters)
        results = [res for res in results if res is not None]

    # Count the number of genes in the filtered clusters
    gene_count = Counter()
    for cluster_id, filtered_genes in results:
        gene_count.update(filtered_genes)

    # Write the results to the output files
    with open(out_file_path, "w") as outfile:
        for cluster_id, filtered_genes in results:
            outfile.write(format_cluster_to_string(cluster_id, filtered_genes))
    write_gene_counts(gene_count, counts_file_path)

    if verbose:
        n_excluded = len(clusters) - len(results)
        print(f"\nPerformed domain filtering on {len(clusters)} tokenised clusters:")
        print(f"  - {len(results)} clusters passed the domain filtering")
        print(
            f"  - {n_excluded} clusters were excluded due to containing less "
            f"than {min_genes} genes with domain hits after domain filtering"
        )
        print(f"The domain-filtered clusters have been saved to {out_file_path}")
