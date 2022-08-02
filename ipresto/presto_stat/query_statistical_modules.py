#!/usr/bin/env python3
"""
Author: Joris Louwen
Student number: 960516530090

Script to query one or more BGCs against a number of statistical modules.
"""

import argparse
from collections import OrderedDict
from functools import partial
from multiprocessing import Pool


def get_commands():
    parser = argparse.ArgumentParser(description="A script to query \
        PRESTO-STAT subclusters to BGCs processed by iPRESTO.")
    parser.add_argument("-i", "--in_file", help="A file that contains BGCs in \
        csv where there is one bgc per line. Each line starts with bgc_name \
        and contains genes represented as domains genes are separated by a \
        ',' and domains in one gene by ';', example: \
        bgc_name1,domain_x;domain_y,domain_a,domain_b", required=True)
    parser.add_argument("-m", "--modfile", dest="modfile", help="Input \
        tsv txt file of modules to query. 6th column should contain \
        modules, header is removed one module per line. Genes are separated \
        by a ',' and domains in one gene by ';'", default=False)
    parser.add_argument('-o', '--out_file', help='Location of output file',
                        required=True)
    parser.add_argument('-c', '--cores', help='Number of cores to use, \
        default = 1', default=1, type=int)
    return parser.parse_args()


def link_mods2bgc(bgc, doms, modules):
    """Returns a tuple of (bgc, [(modules)])

    bgc: string, bgc name
    doms: list of tuples, all domain names from the genes in the bgc
    modules: list of tuples of tuples, each tuple contains the genes of a
        module, each gene contains domains
    """
    modlist = []
    for mod in modules:
        if set(doms).intersection(mod) == set(mod):
            modlist.append(mod)
    return bgc, modlist


def link_all_mods2bgcs(bgcs, modules, cores):
    """Returns a dict of {bgc: [(modules)]}

    bgcs: dict of {bgc: [(domains)])
    modules: list of module tuples of gene tuples
    cores: int, amount of cores to use
    """
    print('\nLinking modules to BGCs')
    pool = Pool(cores, maxtasksperchild=10)
    bgcs_mod = pool.starmap(partial(link_mods2bgc, modules=modules),
                            bgcs.items())
    pool.close()
    pool.join()
    bgc_mod_dict = {pair[0]: pair[1] for pair in bgcs_mod}
    return bgc_mod_dict


def read_clusterfile(infile):
    """Reads a clusterfile into a dictionary of {bgc:[(domains_of_a_gene)]}

    infile: str, filepath
    """
    print("\nReading {}".format(infile))
    with open(infile, 'r') as inf:
        clus_dict = OrderedDict()
        for line in inf:
            line = line.strip().split(',')
            clus = line[0]
            genes = line[1:]
            g_doms = [tuple(gene.split(';')) for gene in genes]
            if clus not in clus_dict.keys():
                clus_dict[clus] = g_doms
            else:
                print("Clusternames not unique, {} read twice".format(clus))
    print("Done. Read {} bgcs".format(len(clus_dict)))
    return clus_dict


def read_mods(modfile):
    """reads modfile to dict of {module:'string_from_file'}

    modfile: filepath

    string from file is stripped
    """
    with open(modfile, 'r') as inf:
        modules = {}
        # {modules:[info]}
        inf.readline()  # header
        for line in inf:
            linesplit = line.strip().split('\t')
            mod = tuple([tuple(gene.split(';')) for gene in
                         linesplit[5].split(',')])
            modules[mod] = line.strip()
        return modules


def write_bgc_mod_fasta(bgcs_with_mods, modinfo, outfile):
    """Writes fasta like file: >BGC\nmod1\nmod2\n

    bgcs_with_mods: {bgc: [(modules)]}
    modinfo: {modules:'string_to_write'}
    outfile: str, filepath
    """
    with open(outfile, 'w') as outf:
        for bgc, moduls in bgcs_with_mods.items():
            outf.write('>{}\n'.format(bgc))
            mod_info_list = [modinfo[mod] for mod in moduls]
            mod_info_list.sort(key=lambda x: int(x[0]))
            for mod_info in mod_info_list:
                outf.write('{}\n'.format(mod_info))


if __name__ == '__main__':
    cmd = get_commands()

    bgcs_main = read_clusterfile(cmd.in_file)
    mods = read_mods(cmd.modfile)
    bgc_with_modlist = link_all_mods2bgcs(bgcs_main, list(mods), cmd.cores)
    print('\nWriting to file')
    write_bgc_mod_fasta(bgc_with_modlist, mods, cmd.out_file)
