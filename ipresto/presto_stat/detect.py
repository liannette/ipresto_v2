from functools import partial
from multiprocessing import Pool
from typing import List
from ipresto.presto_stat.stat_module import StatModule


def map_modules_to_bgc(bgc_id: str, bgc_genes: list, modules: List[StatModule]):
    """
    Returns a tuple of (bgc, [modules]) where each module is a list of genes.

    Parameters:
    bgc_id (str): The name of the biosynthetic gene cluster (BGC).
    bgc_genes (list of tuples): A list of genes, where each gene is a tuple of protein domains.
    modules (list of dicts): A list of dicts where each dict represents a module.

    Returns:
    tuple: A tuple containing the BGC id and a list of module ids that are contained in the BGC.
    """
    contained_modules = []
    for mod in modules:
        # Check if all genes in the module are in the BGC
        if set(mod.tokenised_genes).issubset(set(bgc_genes)):
            contained_modules.append(mod)
    return bgc_id, contained_modules


def detect_modules_in_bgcs(bgcs: dict, modules: List[StatModule], cores: int):
    """Detects STAT subcluster modules in tokenized bgcs.

    Args:
        bgcs (dict): A dictionary where each key is a BGC identifiers and each value
            is a list of genes. Each gene is a tuple of protein domains.
        modules (list): A list of STAT subcluster modules. 
        cores (int): The number of CPU cores to use for parallel processing.

    Returns:
        dict: A dictionary where each key is a BGC identifiers and each value is a list of associated
            modules. Each module is a tuple of genes. Each gene is a tuple of protein domains.
    """
    #todo: change the modules to the dict structure
    pool = Pool(cores, maxtasksperchild=10)
    hits = pool.starmap(
        partial(map_modules_to_bgc, modules=modules), bgcs.items()
    )
    pool.close()
    pool.join()

    return {bgc_id: contained_modules for bgc_id, contained_modules in hits}
