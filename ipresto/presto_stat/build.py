from collections import Counter, defaultdict, OrderedDict
from itertools import product
from multiprocessing import Pool
from functools import partial
from statsmodels.stats.multitest import multipletests
from sympy import binomial as ncr
from math import floor, log10
import networkx as nx

from ipresto.presto_stat.detect import detect_modules_in_bgcs
from ipresto.preprocess.utils import count_non_emtpy_genes
from ipresto.presto_stat.stat_module import StatModule
from typing import List


def remove_infrequent_genes(bgcs, min_genes, min_gene_occurrence, verbose):
    """
    Returns tokenised bgcs with genes replaced  with (-) if they occur < cutoff

    bgcs: dict of {bgc_id:[(domains_in_a_gene)]}
    min_genes: int, minimal non-empty genes (domain combinations) a bgc must have to be included
    verbose: bool, if True print additional info
    cutoff: int, remove genes (domain combinations) that occur less then cutoff
    """
    if verbose:
        print(
            "\nIgnoring infrequent genes that occur less than "
            f"{min_gene_occurrence} times across all BGCs.\nThese genes will not "
            "be considered when generating new PRESTO-STAT subcluster modules."
        )

    # Count occurrences of each gene
    gene_counter = Counter(g for genes in bgcs.values() for g in genes if g != ("-",))
    infreq_genes = {
        g for g, count in gene_counter.items() if count < min_gene_occurrence
    }

    if verbose:
        print(f"{len(infreq_genes)}/{len(gene_counter)} genes are infrequent")

    # Replace infrequent genes with ("-",) and filter out clusters with too few genes
    filtered_bgcs = OrderedDict()
    for bgc_id, genes in bgcs.items():
        filtered_genes = [("-",) if g in infreq_genes else g for g in genes]
        n_filtered_genes = count_non_emtpy_genes(filtered_genes)

        if n_filtered_genes < min_genes:
            if verbose:
                print(
                    f"  Excluding {bgc_id}: only {n_filtered_genes} genes with "
                    f"domain hits after removing infrequent genes (min {min_genes})."
                )
        else:
            filtered_bgcs[bgc_id] = filtered_genes

    if verbose:
        n_removed = len(bgcs) - len(filtered_bgcs)
        print(
            f"{n_removed}/{len(bgcs)} BGCs have been removed due to containing "
            f"less than {min_genes} genes with domain hits after removing "
            "infrequent genes."
        )

    if len(filtered_bgcs) < 2:
        raise ValueError(
            "Not sufficent BGCs remain after filtering out infrequent genes. "
            "Consider lowering the min_gene_occurrence parameter or providing a "
            "larger BGC dataset."
        )

    return filtered_bgcs


def makehash():
    """Function to initialise nested dict"""
    return defaultdict(makehash)


def count_adj(counts, cluster):
    """Counts all adjacency interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    cluster: list of tuples, genes with domains
    """
    if len(cluster) == 1:
        return
    for i, dom in enumerate(cluster):
        if i == 0:
            edge = 1
            adj = [cluster[1]]
        elif i == len(cluster) - 1:
            edge = 1
            adj = [cluster[i - 1]]
        else:
            edge = 2
            adj = [cluster[i - 1], cluster[i + 1]]
            if adj[0] == adj[1] and adj[0] != ("-",):
                # B2 and N2 counts
                prevdom = cluster[i - 1]
                counts[prevdom]["N1"] -= 2
                counts[prevdom]["N2"] += 1
                if dom != ("-",) and dom != prevdom:
                    counts[prevdom]["B1"][dom] -= 2
                    try:
                        counts[prevdom]["B2"][dom] += 1
                    except TypeError:
                        counts[prevdom]["B2"][dom] = 1
        if not dom == ("-",):
            counts[dom]["count"] += 1
            counts[dom]["N1"] += edge
            for ad in adj:
                if ad != ("-",) and ad != dom:
                    try:
                        counts[dom]["B1"][ad] += 1
                    except TypeError:
                        counts[dom]["B1"][ad] = 1


def remove_dupl_doms(cluster):
    """
    Replaces duplicate domains in a cluster with '-', writes domain at the end

    cluster: list of tuples, tuples contain str domain names
    """
    domc = Counter(cluster)
    dupl = [dom for dom in domc if domc[dom] > 1 if not dom == ("-",)]
    if dupl:
        newclus = [("-",) if dom in dupl else dom for dom in cluster]
        for dom in dupl:
            newclus += [("-",), dom]
    else:
        newclus = cluster
    return newclus


def count_coloc(counts, cluster):
    """Counts all colocalisation interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    cluster: list of tuples, genes with domains
    verbose: bool, if True print additional info
    """
    N1 = len(cluster) - 1
    for dom in cluster:
        if not dom == ("-",):
            counts[dom]["count"] += 1
            counts[dom]["N1"] += N1
            coloc = set(cluster)
            try:
                coloc.remove(("-",))
            except KeyError:
                pass
            coloc.remove(dom)
            for colo in coloc:
                try:
                    counts[dom]["B1"][colo] += 1
                except TypeError:
                    counts[dom]["B1"][colo] = 1


def count_interactions(clusdict, verbose):
    """Count all adj and coloc interactions between all domains in clusdict

    clusdict: dict of {cluster:[(gene_with_domains)]}
    verbose: bool, if True print additional info
    Returns two dicts, one dict with adj counts and one with coloc counts
    adj counts:
        { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    coloc counts:
        { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    """
    if verbose:
        print("\nCounting colocalisation and adjacency interactions")
    all_doms = {v for val in clusdict.values() for v in val}
    if ("-",) in all_doms:
        all_doms.remove(("-",))

    # initialising count dicts
    adj_counts = makehash()
    for d in all_doms:
        for v in ["count", "N1", "N2"]:
            adj_counts[d][v] = 0
        for w in ["B1", "B2"]:
            adj_counts[d][w] = makehash()
        # N1: positions adj to one domA, N2: positions adj to two domA
        # B1: amount of domB adj to one domA, B2: positions adj to two domA

    coloc_counts = makehash()
    for d in all_doms:
        for v in ["count", "N1"]:
            coloc_counts[d][v] = 0
        coloc_counts[d]["B1"] = makehash()
        # N1: all possible coloc positions in a cluster, cluster lenght - 1
        # B1: amount of domB coloc with domA

    for clus in clusdict.values():
        count_adj(adj_counts, clus)
        filt_clus = remove_dupl_doms(clus)
        count_coloc(coloc_counts, filt_clus)
    return adj_counts, coloc_counts


def calc_adj_pval(domval_pair, counts, Nall):
    """Returns a list of sorted tuples (domA,domB,pval)

    domval_pair: tuple of (domA, {count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w}} )
    Nall: int, all possible positions
    counts: nested dict { domA:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    """
    domA, vals = domval_pair
    # domains without interactions do not end up in pvals
    if not vals["B1"] and not vals["B2"]:
        return
    pvals = []
    count = vals["count"]
    Ntot = Nall - count
    N1 = vals["N1"]
    N2 = vals["N2"]
    N0 = Ntot - N1 - N2
    interactions = vals["B1"].keys() | vals["B2"].keys()
    for domB in interactions:
        if domB not in vals["B2"]:
            B1 = vals["B1"][domB]
            Btot = counts[domB]["count"]
            pval = float(
                1
                - sum([ncr(N0, (Btot - d)) * ncr(N1, d) for d in range(B1)])
                / ncr(Ntot, Btot)
            )
        elif vals["B1"][domB] == 0:
            B2 = vals["B2"][domB]
            Btot = counts[domB]["count"]
            pval = float(
                1
                - sum([ncr(N0, (Btot - d)) * ncr(N2, d) for d in range(B2)])
                / ncr(Ntot, Btot)
            )
        else:
            B1 = vals["B1"][domB]
            B2 = vals["B2"][domB]
            Btot = counts[domB]["count"]
            pval = float(
                1
                - sum(
                    [
                        ncr(N0, Btot - d1 - d2) * ncr(N1, d1) * ncr(N2, d2)
                        for d1, d2 in product(range(B1 + 1), range(B2 + 1))
                        if d1 + d2 != B1 + B2
                    ]
                )
                / ncr(Ntot, Btot)
            )
        ab_int = sorted((domA, domB))
        pvals.append((ab_int[0], ab_int[1], pval))
    return pvals


def calc_adj_pval_wrapper(count_dict, clusdict, cores, verbose):
    """Returns list of tuples of corrected pvals for each gene pair

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    clusdict: dict of {cluster:[(domains_in_a_gene)]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    """
    if verbose:
        print("Calculating adjacency p-values")
    N = sum([len(values) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=5)
    pvals_ori = pool.map(
        partial(calc_adj_pval, counts=count_dict, Nall=N), count_dict.items()
    )
    pool.close()
    pool.join()
    # remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    # to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0], tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0], p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0], p[1])] != 2]
            print("  error with domain pairs {}".format(", ".join(p_excl)))
            print("  these are excluded")
    # Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method="fdr_by")[1]
    # adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in zip(
        pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]
    ):
        assert ab1[0] == ab2[0] and ab1[1] == ab2[1]
        pmax = max([p1, p2])
        ptups.append(((ab1[0], ab1[1]), pmax))
    return ptups


def calc_coloc_pval(domval_pair, counts, Nall):
    """Returns a list of sorted tuples (domA,domB,pval)

    domval_pair: tuple of (domA, { count:x,N1:y,B1:{dom2:v,dom3:w } })
    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    Nall: int, all possible positions in all clusters
    """
    domA, vals = domval_pair
    # domains without interactions do not end up in pvals
    if not vals["B1"]:
        return
    pvals = []
    count = vals["count"]
    Ntot = Nall - count
    N1 = vals["N1"]
    N0 = Ntot - N1
    interactions = vals["B1"].keys()
    for domB in interactions:
        B1 = vals["B1"][domB]
        Btot = counts[domB]["count"]
        pval = float(
            1
            - sum([ncr(N0, (Btot - d)) * ncr(N1, d) for d in range(B1)])
            / ncr(Ntot, Btot)
        )
        ab_int = sorted((domA, domB))
        pvals.append((ab_int[0], ab_int[1], pval))
    return pvals


def calc_coloc_pval_wrapper(count_dict, clusdict, cores, verbose):
    """Returns list of tuples of corrected pvals for each domain pair

    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    clusdict: dict of {cluster:[domains]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    """
    print("Calculating colocalisation p-values")
    N = sum([len(remove_dupl_doms(values)) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=1)
    pvals_ori = pool.map(
        partial(calc_coloc_pval, counts=count_dict, Nall=N), count_dict.items()
    )
    pool.close()
    pool.join()
    # remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    # to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0], tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0], p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0], p[1])] != 2]
            print("  error with domain pairs {}".format(", ".join(p_excl)))
            print("  these are excluded")
    # Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method="fdr_by")[1]
    # adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in zip(
        pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]
    ):
        assert ab1[0] == ab2[0] and ab1[1] == ab2[1]
        pmax = max([p1, p2])
        ptups.append(((ab1[0], ab1[1]), pmax))
    return ptups


def keep_lowest_pval(colocs, adjs):
    """
    Returns all domain pairs with their lowest pvalue as an edge for nx

    colocs, adjs: list of tuples [((dom1,dom2),pval)]
    Tuples look like (dom1,dom2,{pval:x})
    """
    pvals = colocs + adjs
    counter = Counter(list(zip(*pvals))[0])
    dupl = sorted([tup for tup in pvals if counter[tup[0]] == 2])
    uniques = [tup for tup in pvals if counter[tup[0]] == 1]
    lowest = []
    for p1, p2 in zip(dupl[::2], dupl[1::2]):
        pmin = min([p1[1], p2[1]])
        lowest.append((p1[0][0], p1[0][1], {"pval": pmin}))
    uniques = [(tup[0][0], tup[0][1], {"pval": tup[1]}) for tup in uniques]
    return lowest + uniques


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


def find_modules_by_pval_cutoff(pval_cutoff, edges):
    """
    Returns modules found given a specific pval cutoff as (pval_cutoff, {modules})

    pval_cutoff: float, cutoff for detecting modules
    gene_pairs: list of tuples, ('gene1', 'gene2', {"pval": pvalue})
    Modules are all maximal cliques with length > 2
    """
    filtered_edges = (e for e in edges if e[2]["pval"] <= pval_cutoff)
    module_graph = generate_graph(filtered_edges, False)
    cliques = nx.algorithms.clique.find_cliques(module_graph)
    modules = {tuple(sorted(module)) for module in cliques if len(module) > 2}
    return pval_cutoff, modules


def round_to_n(x, n):
    """Round x to n significant decimals

    x: int/float
    n: int
    """
    if x <= 0:
        return 0
    return round(x, -int(floor(log10(x))) + (n - 1))


def identify_significant_modules(edges, max_pval, cores, verbose):
    """
    Generates a dictionary of modules with their strictest p-value cutoff.

    This function identifies all modules with a p-value lower than the specified
    cutoff and returns a dictionary where the keys are the modules and the values
    are the strictest p-value cutoff for each module.

    Args:
        edges (list of tuples): A list of tuples containing two genes and their
            associated interaction p-values: ('gene1', 'gene2', {"pval": pvalue}))
        max_pval (float): The p-value cutoff for significance.
        cores (int): The number of CPU cores to use for parallel processing.
        verbose (bool): If True, prints additional information during execution.

    Returns:
        list: A list of dictionaries, with each dictionary representing a module
            - "module_id" (int): A unique identifier for the module.
            - "n_genes" (int): The number of genes in the module.
            - "n_domains" (int): The total number of domains in the module.
            - "strictest_pval" (float): The strictest p-value cutoff for the module.
            - "module" (tuple): The module itself.
    """
    if verbose:
        print(
            f"\nIdentifying subcluster modules applying a maximum interaction "
            f"p-value cutoff of {max_pval}"
        )

    significant_edges = [e for e in edges if e[2]["pval"] <= max_pval]
    if verbose:
        print(f"  {len(significant_edges)} significant gene pair interactions")

    pval_cutoffs = {pv["pval"] for pv in list(zip(*significant_edges))[2]}
    if len(pval_cutoffs) > 100000:  # reduce the number of pvals to loop through
        pval_cutoffs = {round_to_n(x, 3) for x in pval_cutoffs}
    if verbose:
        print(f"  looping through {len(pval_cutoffs)} p-value cutoffs")

    pool = Pool(cores, maxtasksperchild=10)
    results = pool.imap(
        partial(find_modules_by_pval_cutoff, edges=significant_edges),
        pval_cutoffs,
        chunksize=250,
    )
    pool.close()
    pool.join()

    modules = {}
    module_id = 1
    for result in results:
        pval_cutoff, generated_modules = result
        for mod in generated_modules:
            if mod not in modules:
                modules[mod] = StatModule(mod, pval_cutoff, module_id)
                module_id += 1
            else:
                strictest_pval = min(modules[mod].strictest_pval, pval_cutoff)
                modules[mod].strictest_pval = strictest_pval

    print(f"Detected {len(modules)} modules.")
    return list(modules.values())


def calculate_interaction_pvals(
    bgcs,
    cores,
    verbose,
):
    """
    Returns:
        p_values (list of tuples): A list of tuples containing two genes and their
                associated interaction p-values: ('gene1', 'gene2', {"pval": pvalue}))
    """
    adj_counts, c_counts = count_interactions(bgcs, verbose)
    adj_pvals = calc_adj_pval_wrapper(adj_counts, bgcs, cores, verbose)
    col_pvals = calc_coloc_pval_wrapper(c_counts, bgcs, cores, verbose)
    pvals = keep_lowest_pval(col_pvals, adj_pvals)
    # todo: keep from crashing when there are no significant modules
    return pvals


def count_module_occurrences(modules, bgcs, cores):
    """Counts the number of occurences of each module in the BGC dataset.

    Args:
        modules (list): A list of dictionaries representing modules.
        bgcs (dict): A dictionary where each key is a BGC identifier and each value
            is a list of genes. Each gene is a tuple of protein domains.
        cores (int): The number of CPU cores to use for parallel processing.

    Returns:
        dict: A dictionary where each key is a module identifier and each value
            is the number of occurrences of that module in the BGC dataset.
    """
    detected_modules = detect_modules_in_bgcs(bgcs, modules, cores)
    module_counts = Counter(
        mod for bgc, mods in detected_modules.items() for mod in mods
    )
    for mod in modules:
        mod.occurences = module_counts.get(mod, 0)

    return modules


def generate_stat_modules(
    out_dir,
    bgcs,
    min_genes,
    max_pval,
    min_gene_occurence,
    cores,
    verbose,
):
    # remove genes that occur too infrequent in the bgc dataset
    bgcs = remove_infrequent_genes(bgcs, min_genes, min_gene_occurence, verbose)

    # find the modules
    p_values = calculate_interaction_pvals(bgcs, cores, verbose)
    modules = identify_significant_modules(p_values, max_pval, cores, verbose)

    # count the number of occurences of each module in the BGC dataset
    modules = count_module_occurrences(modules, bgcs, cores)

    return modules


def filter_infrequent_modules(
    modules: List[StatModule], min_occurence: int
) -> List[StatModule]:
    """Removes modules that occur less than twice in the BGC dataset.

    Args:
        modules (list): A list of dictionaries representing modules.

    Returns:
        list: A filtered list of dictionaries representing modules that occur more than once.
    """
    filtered_modules = [mod for mod in modules if mod.occurences > 1]
    for new_id, module in enumerate(filtered_modules, start=1):
        module.module_id = new_id
    return filtered_modules
