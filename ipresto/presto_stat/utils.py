import csv
from collections import OrderedDict
from ipresto.preprocess.utils import count_non_emtpy_genes
from ipresto.presto_stat.stat_module import StatModule
from typing import List


def read_clusterfile(clusterfile, m_gens, verbose):
    """
    Reads a cluster file into a dictionary of {bgc: [(domains_of_a_gene)]}.

    Parameters:
        clusterfile (str): The file path to the input cluster file.
        m_gens (int): The minimum number of genes with domains a cluster should have.
        verbose (bool): If True, print additional information during processing.

    Returns:
        dict: A dictionary where keys are cluster names (bgc) and values are lists of tuples,
            each tuple representing the protein domains in a gene.

    The function also prints the number of clusters read and the number of clusters excluded
    due to having fewer genes than the specified minimum. Clusters with non-unique names
    will trigger a warning message.

    Notes:
    - Clusters with fewer than `m_gens` genes are not included in the returned dictionary.
    - The function also prints the number of clusters excluded due to having fewer genes than `m_gens`.
    - Domains represented by '-' are not counted in the gene count.
    """
    print(f"\nReading {clusterfile}")
    with open(clusterfile, "r") as inf:
        clusters = OrderedDict()
        for line in inf:
            line = line.strip().split(",")
            clus = line[0]
            genes = line[1:]
            g_doms = [tuple(gene.split(";")) for gene in genes]
            if clus not in clusters.keys():
                clusters[clus] = g_doms
            else:
                print(f"Warning: Duplicate cluster ID {clus}.")

    # remove clusters with too few genes
    if m_gens > 0:
        clusters = remove_empty_clusters(clusters, m_gens, verbose)

    print(f"Done. Read {len(clusters)} clusters")
    return clusters


def remove_empty_clusters(clusters, min_genes, verbose):
    """
    Removes clusters with fewer than `min_genes` genes from the clusters dictionary.

    Parameters:
        clusters (dict): The dictionary of clusters to filter.
        min_genes (int): The minimum number of genes with domains a cluster should have.
        verbose (bool): If True, print additional information during processing.

    Returns:
        dict: The filtered dictionary of clusters, where each cluster has at least `min_genes` genes.

    """
    to_delete = []
    for clus, genes in clusters.items():
        n_genes = count_non_emtpy_genes(genes)
        if n_genes < min_genes:
            if verbose:
                print(
                    f"  Excluding {clus}: contains only {n_genes} genes with "
                    f"domain hits."
                )
            to_delete.append(clus)

    for clus in to_delete:
        del clusters[clus]

    if to_delete:
        print(
            f"{len(to_delete)} clusters excluded for having less than "
            f"{min_genes} genes with domain hits (minimum required: {min_genes})."
        )

    return clusters


def tokenized_genes_to_string(tokenized_genes):
    """
    Converts a list of tokenized genes (tuples of tuples) into a string representation.

    Parameters:
        tokenized_genes (list of tuples): A list of tokenized genes, where each gene is
            represented as a tuple of tuples.

    Returns:
        str: A string representation of the tokenized genes, with each gene separated by a comma.
    """
    return ",".join([";".join(gene) for gene in tokenized_genes])


def string_to_tokenized_genes(genes_string):
    """
    Converts a string representation of tokenized genes back into a list of tuples.

    Parameters:
        genes_string (str): A string representation of tokenized genes, where each gene is
            separated by a comma and each tuple is separated by a semicolon.

    Returns:
        list of tuples: A list of tokenized genes, where each gene is represented as a tuple of tuples.
    """
    return [tuple(gene.split(";")) for gene in genes_string.split(",")]


def write_stat_modules(modules: List[StatModule], file_path: str):
    """
    Writes the statistical modules to a file in tab-separated format.

    Parameters:
        modules (List[dict]): A list of dictionaries where each dictionary represents a module.
        file_path (str): The file path to the output file.

    The function writes the header based on the keys of the first dictionary in the list,
    and then writes each module's values in tab-separated format.
    """
    if not modules:
        print("No modules to write.")
        return

    header = modules[0].to_dict().keys()
    with open(file_path, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for module in modules:
            row = module.to_dict()
            row["tokenised_genes"] = tokenized_genes_to_string(row["tokenised_genes"])
            writer.writerow(row)


def read_stat_modules(file_path):
    """
    Reads statistical modules from a tab-separated file into a list of dictionaries.

    Parameters:
        file_path (str): The file path to the input file.

    Returns:
        list of dict: A list of dictionaries where each dictionary represents a module.
    """
    with open(file_path, "r", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        modules = []
        for row in reader:
            tokenised_genes = string_to_tokenized_genes(row["tokenised_genes"])
            strictest_pval = float(row["strictest_pval"])
            modules.append(
                StatModule(tokenised_genes, strictest_pval, row["module_id"])
            )
        return modules


def write_detected_stat_module_ids(detected_modules, out_file_path):
    with open(out_file_path, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        header = ["bgc_id", "module_ids"]
        writer.writerow(header)  # Write the header
        for bgc_id, modules in detected_modules.items():
            module_ids = [str(mod.module_id) for mod in modules]
            writer.writerow([bgc_id, ",".join(module_ids)])


def write_detected_stat_modules(detected_modules, out_file_path):
    with open(out_file_path, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for bgc_id, modules in detected_modules.items():
            outfile.write(f">{bgc_id}\n")
            for module in modules:
                row = module.to_dict()
                tokenised_genes = tokenized_genes_to_string(row["tokenised_genes"])
                writer.writerow(
                    [
                        row["module_id"],
                        tokenised_genes,
                    ]
                )
