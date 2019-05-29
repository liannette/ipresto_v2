#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert the genes in BGCs into strings of domains, filter BGCs
based on similarity, and detect sub-clusters based on a statistical method
and an LDA alorithm.

Usage:
python3 gene_module_detection_pipeline.py -h

Example usage:
python3 gene_module_detection_pipeline.py -i ./testdata -o ./testdata_domains
--hmm_path ./domains/Pfam_100subs_tc.hmm --exclude final -c 12 -e True

Notes:
Only handles gbk files with one cluster

Layout:
get_commands
process_gbks
convert_gbk2fasta
run_hmmscan
hmmscan_wrapper
parse_domtab
sign_overlap
parse_dom_wrapper
read_clusterfile
calc_adj_index
is_contained
generate_edges
generate_edge
generate_graph
find_representatives
find_all_representatives
write_filtered_bgcs
remove_infr_doms
remove_dupl_doms
count_adj
count_coloc
makehash
count_interactions
calc_adj_pval_wrapper
calc_adj_pval
calc_coloc_pval_wrapper
keep_lowest_pval
visualise_graph
generate_modules_wrapper
generate_modules
write_module_file

Required:
python 3
hmmscan
'''

import argparse
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter, defaultdict
from copy import deepcopy
from functools import partial
from glob import glob, iglob
from itertools import combinations, product
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import networkx as nx
from operator import itemgetter
import os
import random
import re
from statsmodels.stats.multitest import multipletests
import subprocess
from sympy import binomial as ncr
import time

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn bgcs in \
        gbk files into strings of domains using a domain hmm database and to \
        reduce redundancy by filtering out similar bgcs.")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        directory of gbk files", required=True)
    parser.add_argument("--exclude", dest="exclude", default=["final"],
        nargs="+", help="If any string in this list occurs in the gbk \
        filename, this file will not be used for the analysis. \
        (default: final)")
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("--hmm_path", dest="hmm_path", required=True,
        help="File containing domain hmms that is hmmpress-processed.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    parser.add_argument("-v", "--verbose", dest="verbose", required=False,
        action="store_true", default=False, help="Prints more detailed \
        information.")
    parser.add_argument("-d", "--domain_overlap_cutoff", 
        dest="domain_overlap_cutoff", default=0.1, help="Specify at which \
        overlap percentage domains are considered to overlap. Domain with \
        the best score is kept (default=0.1).")
    parser.add_argument("-e", "--exclude_contig_edge",
        dest="exclude_contig_edge", default=True, type=bool, help="\
        Exclude clusters that lie on a contig edge")
    parser.add_argument("-m", "--min_genes", dest="min_genes", default=0,
        help="Provide the minimum size of a BGC to be included in the \
        analysis. Default is 0 genes", type=int)
    parser.add_argument("--min_doms", dest="min_doms", default=0,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 0 domains", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    parser.add_argument("-p", "--pval_cutoff", dest="pval_cutoff", \
        default = 0.1, type=float, help='P-value cutoff for determining a \
        significant interaction in module detection (default: 0.1)')
    parser.add_argument("--use_fastas", dest="use_fastas", default=None, \
        help="Use already created fasta files from some folder")
    parser.add_argument("--use_domtabs", dest="use_domtabs", default=None, \
        help="Use already created domtables from some folder")
    parser.add_argument("--include_list", dest="include_list", default=None, \
        help="If provided only the domains in this file will be taken into \
        account in the analysis. One line should contain one Pfam ID \
        (default: None - meaning all Pfams from database)")
    parser.add_argument("--start_from_clusterfile", default=None, help="A file\
        with BGCs and domain-combinations to start with (csv and domains in a\
        gene separated by ';'). This overwrites in_folder (which still has to\
        be supplied symbolically) and use_domtabs/use_fastas.")
    return parser.parse_args()

def process_gbks(input_folder, output_folder, exclude, exclude_contig_edge,\
    min_genes, cores, verbose, existing_fasta_folder):
    '''Convert gbk files from input folder to fasta files for each gbk file

    input_folder, outpu_folder: str
    exclude: list of str, files will be excluded if part of the file name
        is present in this list
    exclude_contig_edge: bool
    min_genes: int
    verbose: bool, print additional info to stdout
    '''
    if not os.path.isdir(output_folder):
        subprocess.check_call("mkdir {}".format(output_folder), shell = True)
    if input_folder.endswith('/'):
        base, inf = os.path.split(input_folder[:-1])
    else:
        base, inf = os.path.split(input_folder)
    out_fasta = os.path.join(output_folder, inf+'_fasta')
    if not os.path.isdir(out_fasta):
        subprocess.check_call("mkdir {}".format(out_fasta), shell = True)
    print("\nProcessing gbk files into fasta files.")
    files = iglob(os.path.join(input_folder, "*.gbk"))
    done = []
    pool = Pool(cores, maxtasksperchild=20)
    for file_path in files:
        pool.apply_async(convert_gbk2fasta, args=(file_path, out_fasta,\
            exclude_contig_edge, min_genes, exclude, verbose, \
            existing_fasta_folder), callback=lambda x: done.append(x))
    pool.close()
    pool.join()
    processed = [val for val in done if val]
    fastas_in_existing_folder = [val for val in processed if type(val)==str]
    excluded = len([val for val in done if val == False])
    filtered = len([val for val in done if val == None])
    print("Processed {} gbk files into {} fasta files.".format(\
        len(processed)+excluded+filtered, len(processed)))
    print(" excluded {} files containing {}".format(excluded,\
        ' or '.join(exclude)))
    print(" filtered {} files that didn't pass constraints".format(\
        filtered))
    return out_fasta, fastas_in_existing_folder

def convert_gbk2fasta(file_path, out_folder, exclude_contig_edge, min_genes,\
    exclude, verbose, existing_fasta_folder):
    '''Convert one gbk file to a fasta file in out_folder

    file_path, out_folder: strings
    exclude_contig_edge: bool
    min_genes: int
    verbose: bool, print additional info to stdout

    Returns True for a successful conversion to fasta, None if there is a
    contig edge or min_genes is not passed. False is returned if any of
    the exclude list is in the filename

    If the fasta file already exists in the existing_fasta_folder, the fasta
    file is not created, True is returned
    '''
    file_name = os.path.split(file_path)[1]
    if any([word in file_name for word in exclude]):
        return False
    name = file_name.split('.gbk')[0]
    name_extend = '{}.fasta'.format(name)
    outfile = os.path.join(out_folder, name_extend)
    if existing_fasta_folder:
        existing_file = os.path.join(existing_fasta_folder, name_extend)
    else:
        existing_file = 'not_an_existing_file'
    seqs = OrderedDict()
    num_genes = 0
    if not os.path.isfile(outfile) and not os.path.isfile(existing_file):
        try:
            record = next(SeqIO.parse(file_path, 'genbank'))
        except ValueError as e:
            print(" Excluding {}: {}".format(file_path, e))
            return
        for feature in record.features:
            if feature.type == 'cluster':
                if "contig_edge" in feature.qualifiers:
                    if feature.qualifiers["contig_edge"][0] == "True":
                        if exclude_contig_edge:
                            if verbose:
                                print("  excluding {}: {}".format(file_name,\
                                    "contig edge"))
                            return
            if feature.type == 'CDS':
                gene_id = ""
                if "gene" in CDS.qualifiers:
                    gene_id = CDS.qualifiers.get('gene',"")[0]
                protein_id = ""
                if "protein_id" in CDS.qualifiers:
                    protein_id = CDS.qualifiers.get('protein_id',"")[0]
                head = '_'.join(name,gene_id,protein_id)
                header = ">{}_{}".format(head, num_genes+1)
                seqs[header] = feature.qualifiers['translation'][0]
                if seqs[header] == '':
                    print('  {} does not have a translation'.format(header))
                num_genes +=1
        if num_genes < min_genes:
            if verbose:
                print("  excluding {}: less than {} genes".format(file_path,\
                    min_genes))
            return
        with open(outfile, 'w') as out:
            for seq in seqs:
                compl_header = '{}/{}'.format(seq,num_genes)
                out.write('{}\n{}\n'.format(compl_header, seqs[seq]))
    elif not os.path.exists(outfile):
            return existing_file
    return True

def run_hmmscan(fasta_file, hmm_file, out_folder, verbose, \
    existing_dom_folder):
    """
    Runs hmmscan on fasta file to generate a domtable file

    fasta_file, hmm_file, out_folder: strings of file paths
    verbose: bool
    """
    if os.path.isfile(fasta_file):
        name = os.path.split(fasta_file)[1].split('.fasta')[0]
        new_name = name+".domtable"
        out_name = os.path.join(out_folder, new_name)
        if existing_dom_folder:
            existing_out_name = os.path.join(existing_dom_folder, new_name)
        else:
            existing_out_name = './not_an_existing_file'
        log = os.path.join(out_folder, 'hmmlog.txt')
        if not os.path.isfile(out_name) and not \
            os.path.isfile(existing_out_name):
            hmmscan_cmd = (\
                "hmmscan -o {} --cpu 0 --domtblout {} --cut_tc {} {}".format(\
                log, out_name, hmm_file, fasta_file))
            if verbose:
                print("  " + hmmscan_cmd)
            subprocess.check_call(hmmscan_cmd, shell=True)
        elif os.path.isfile(out_name):
            if verbose:
                print("  {} existed. hmmscan not run again".format(out_name))
        else:
            #return the existing file
            if verbose:
                print("  {} existed. hmmscan not run again".format(\
                    existing_out_name))
            return existing_out_name
    else:
        raise SystemExit("Error running hmmscan: {} doesn't exist".format(\
            fasta_file))
    return False

def hmmscan_wrapper(input_folder, hmm_file, verbose, cores, \
    existing_fasta_files, existing_dom_folder):
    '''
    Runs hmmscan on all fasta files in input_folder with hmm_file as hmm db

    input_folder, hmm_file: strings of file paths
    verbose: bool, if True print additional information
    cores: int, amount of cores to use
    existing_fasta_files: list of str, filepaths of existing fastas to use
    existing_dom_folder: str, file path to a folder to see if a domtab
        already exists
    '''
    if input_folder.endswith('/'):
        out_folder = input_folder[:-7]+'_domtables'
    else:
        out_folder = input_folder[:-6]+'_domtables'
    if not os.path.isdir(out_folder):
        subprocess.check_call("mkdir {}".format(out_folder), shell = True)
    print("\nRunning hmmscan on fastas to generate domtables.")
    if existing_fasta_files:
        files = glob(os.path.join(input_folder, "*.fasta")) + \
            existing_fasta_files
    else:
        files = glob(os.path.join(input_folder, "*.fasta"))
    pool = Pool(cores, maxtasksperchild=20)
    #maxtasksperchild=1:each process respawns after completing 1 chunk
    done = []
    for file_path in files:
        pool.apply_async(run_hmmscan,args=(file_path, hmm_file,
            out_folder, verbose, existing_dom_folder),
            callback=lambda x: done.append(x))
    pool.close()
    pool.join() #make the code in main wait for the pool processes to finish
    domtabs_in_existing_folder = [val for val in done if val]
    print("Processed {} fasta files into domtables.".format(\
        len(files)))
    return out_folder, domtabs_in_existing_folder

def parse_domtab(domfile, clus_file, sum_file, min_overlap, verbose):
    '''Parses domtab into a cluster domain file (csv)

    domfile: string, file path
    clus_file: open file for writing
    sum_file: open file for writing
    min_overlap : float, the amount of overlap two domains must have for it
        to be considered overlap

    clus_file will look like this:
    Clus1,dom1;dom2,-(gene without domain)\\nClus2,dom1..
    Genes are separated by commas while domains in the genes are separated
        by ;
    '''
    if verbose:
        print("  parsing domtable {}".format(domfile))
    queries = SearchIO.parse(domfile, 'hmmscan3-domtab')
    cds_before = 0
    #list of lists for the domains in the cluster where each sublist is a gene
    cluster_doms = []
    #for every cds that has a hit
    for query in queries:
        dom_matches = []
        cds_num, total_genes = map(int,query.id.split('_')[-1].split('/'))
        #for every hit in each cds
        for hit in query:
            match = hit[0]
            domain = match.hit_id
            range_q = match.query_range
            bitsc = match.bitscore
            dom_matches.append((domain, range_q, bitsc))
        dels = []
        if len(query) > 1:
            for i in range(len(query)-1):
                for j in range(i+1, len(query)):
                    #if there is a significant overlap delete the one with
                    #the lower bitscore
                    if sign_overlap(dom_matches[i][1],dom_matches[j][1],
                        min_overlap):
                        if dom_matches[i][2] >= dom_matches[j][2]:
                            dels.append(j)
                        else:
                            dels.append(i)
        cds_doms = tuple(dom_matches[i][0] for i in range(len(query)) \
            if i not in dels)
        #if a cds has no domains print '-' in output
        gene_gap = cds_num - cds_before -1
        if gene_gap > 0:
            gaps = [('-',) for x in range(gene_gap)]
            cluster_doms += gaps
        cluster_doms.append(cds_doms)
        cds_before = cds_num
    if cds_before == 0:
        print(' excluding {} no domain hits present'.format(domfile))
        return
    end_gap = total_genes - cds_num
    if end_gap > 0:
        gaps = [('-',) for x in range(end_gap)]
        cluster_doms += gaps
    clus_file.write('{},{}\n'.format(\
        os.path.split(domfile)[-1].split('.domtable')[0],
        ','.join(';'.join(gene) for gene in cluster_doms)))
    return cluster_doms

def sign_overlap(tup1, tup2, cutoff):
    '''
    Returns true if there is an overlap between two ranges higher than cutoff

    tup1, tup2: tuples of two ints, start and end of alignment
    cutoff: float, fraction that two alignments are allowed to overlap

    Overlap is be calculated with the smallest domain alignment to be strict
    '''
    overlap = len(range(max(tup1[0], tup2[0]), min(tup1[1], tup2[1])))
    if overlap > 0:
        if overlap > min(abs(tup1[0]-tup1[1]), abs(tup2[0]-tup2[1]))*cutoff:
            return True
    return False

def parse_dom_wrapper(in_folder, out_folder, cutoff, verbose, \
    domtabs_in_existing_folder):
    '''Calls parse_domtab on all domtable files to create a clusterfile

    in_folder, out_folder: strings, filepaths
    cutoff: float, cutoff value for domain overlap
    '''
    print("\nParsing domtables from folder {}".format(in_folder))
    if domtabs_in_existing_folder:
        domtables = glob(os.path.join(in_folder, '*.domtable')) + \
            domtabs_in_existing_folder
    else:
        domtables = iglob(os.path.join(in_folder, '*.domtable'))
    in_name = os.path.split(in_folder)[1].split('_domtables')[0]
    out_file = os.path.join(out_folder, in_name+'_clusterfile.csv')
    sumfile = os.path.join(out_folder, in_name+'_dom_hits.txt')
    stat_file = os.path.join(out_folder, in_name+'_domstats.txt')
    domc = Counter()
    with open(out_file, 'w') as out, open(sumfile,'w') as sumf:
        for domtable in domtables:
            doms = parse_domtab(domtable, out, sumf, cutoff, verbose)
            domc.update(doms)
    with open(stat_file, 'w') as stat:
        stat.write("#Total\t{}\n".format(sum(domc.values())))
        for dom, count in domc.most_common():
            stat.write("{}\t{}\n".format(';'.join(dom),count))
    print("Result in {}".format(out_file))
    print(" statistics about doms in {}".format(stat_file))
    return out_file

def read_clusterfile(infile, m_gens, verbose):
    """Reads a clusterfile into a dictionary of {bgc:[(domains_of_a_gene)]}

    infile: str, filepath
    m_gens: int, minimum of genes with domains a cluster should have
    verbose: bool, if True print additional info

    clusters with less than m_gens genes are not returned
    It also returns a dict {bgc:amount_of_domains} where -'s are not counted
    and all domains of all genes are added to an int
    """
    print("\nReading {}".format(infile))
    filtered = 0
    with open(infile, 'r') as inf:
        clus_dict = OrderedDict()
        for line in inf:
            line = line.strip().split(',')
            clus = line[0]
            genes = line[1:]
            g_doms = [tuple(gene.split(';')) for gene in genes]
            if len([g for g in genes if g != ('-',)]) < m_gens:
                filtered +=1
                if verbose:
                    print("  excluding {} less than min genes".format(clus))
                continue
            if not clus in clus_dict.keys():
                clus_dict[clus] = g_doms
            else:
                print("Clusternames not unique, {} read twice".format(clus))
    print("Done. Read {} clusters".format(len(clus_dict)))
    print(" {} clusters have less than {} genes and are excluded".format(\
        filtered,m_gens))
    return clus_dict

def calc_adj_index(clus1, clus2):
    '''Returns the adjacency index between two clusters

    clus1, clus2: list of str, domainlist of a cluster

    If there is an empty gene between two domains these two domains are not
        adjacent
    '''
    #generate all unique domain pairs
    dom_p1 = {tuple(sorted(dp)) for dp in zip(*(clus1[:-1],clus1[1:])) \
        if not '-' in dp}
    dom_p2 = {tuple(sorted(dp)) for dp in zip(*(clus2[:-1],clus2[1:])) \
        if not '-' in dp}
    #if doms are separated by '-' then there are no dom pairs. if happens ai=0
    if not dom_p1 or not dom_p2:
        return 0.0        
    ai = len(dom_p1 & dom_p2)/len(dom_p1 | dom_p2)
    return ai

def is_contained(clus1, clus2):
    '''
    Returns a bool if all domains from one of the clusters are in the other

    clus1, clus2: list of str, domainlist of a cluster
    '''
    one_in_two = all(dom in clus2 for dom in clus1 if not dom == '-')
    two_in_one = all(dom in clus1 for dom in clus2 if not dom == '-')
    if one_in_two or two_in_one:
        return True
    return False

def generate_edges(dom_dict, cutoff, cores):
    '''Returns a pair of clusters in a tuple if ai/contained above cutoff

    dom_dict: dict {clus1:[domains]}, clusters linked to domains
    cutoff: float, between 0-1, when clusters are similar
    cores: int, amount of cores used for calculation
    '''
    print("\nGenerating similarity scores")
    clus_names = list(dom_dict.keys())
    pairs = combinations(clus_names, 2)
    pool = Pool(cores, maxtasksperchild = 100)
    edges = pool.map(partial(generate_edge, d_dict = dom_dict, \
        cutoff = cutoff), pairs)
    edges = [edge for edge in edges if not edge == None]
    print("Done. {} pairs above threshold".format(len(edges)))
    return edges

def generate_edge(pair, d_dict, cutoff):
    '''
    Calculate similarity scores between two bgcs and return if above cutoff

    pair: tuple of 2 strings, 2 clusternames
    d_dict: dict of {clustername:domains}
    cutoff: float
    A tuple is returned that can be read as an edge by nx.Graph.add_edges_from
    '''
    p1,p2 = pair
    #unpack the domains from the gene tuples
    clus1 = [d for dom in d_dict[p1] for d in dom]
    clus2 = [d for dom in d_dict[p2] for d in dom]
    contained = is_contained(clus1, clus2)
    ai = calc_adj_index(clus1, clus2)
    if ai == None:
        print('  error in generate_edge: {}'.forat(pair))
    if contained or ai > cutoff:
        # print(pair,ai,contained)
        return(p1,p2,{'ai':ai,'contained':contained})

def generate_graph(edges, verbose):
    '''Returns a networkx graph

    edges: list of tuples, (pair1,pair2,{attributes})
    '''
    g = nx.Graph()
    g.add_edges_from(edges)
    if verbose:
        print('\nGenerated graph with:')
        print(' {} nodes'.format(g.number_of_nodes()))
        print(' {} edges'.format(g.number_of_edges()))
    return g

def find_representatives(clqs, d_l_dict, graph):
    '''
    Returns {representative:[clique]} based on bgc with most domains in clique

    clqs: list of lists of strings, cliques of clusters
    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    graph: networkx graph structure of the cliques
    The longest cluster is chosen (most domains). If there are multiple
        longest clusters then the cluster with the least connections is
        chosen (to preserve most information).
    '''
    reps_dict = OrderedDict()
    dels = set() #set of nodes for which a representative has been found
    for cliq in clqs:
        cliq = [clus for clus in cliq if not clus in dels]
        if cliq:
            domlist = [(clus,d_l_dict[clus]) for clus in cliq]
            maxdoml = max(doms[1] for doms in domlist)
            clus_maxlen = [clus for clus, doml in domlist \
                if doml == maxdoml]
            if len(clus_maxlen) > 1:
                min_degr = min([deg for clus, deg in \
                    graph.degree(clus_maxlen)])
                rep = random.choice([clus for clus in clus_maxlen \
                    if graph.degree(clus) == min_degr])
            else:
                rep = clus_maxlen[0]
            try:
                reps_dict[rep].update(cliq)
            except KeyError:
                reps_dict[rep] = set(cliq)
            cliq.remove(rep)
            dels.update(cliq)
    return reps_dict

def find_all_representatives(d_l_dict, g):
    '''Iterates find_representatives until there are no similar bgcs

    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    g: networkx graph structure containing the cliques
    all_reps_dict: dict of {representative:[represented]}
    '''
    print('\nFiltering out similar bgcs.')
    all_reps_dict = {}
    subg = g.subgraph(g.nodes)
    i = 0
    while subg.number_of_edges() != 0:
        print(\
        '  iteration {}, edges (similarities between bgcs) left: {}'.format(\
            i,subg.number_of_edges()))
        cliqs = nx.algorithms.clique.find_cliques(subg)
        #make reproducible by making the cliqs have the same order every time
        #sort first each cliq alphabetically, then cliqs alphabetically,
        #then on length, so longest are first and order is the same
        cliqs = sorted(sorted(cl) for cl in cliqs if len(cl) > 1)
        cliqs.sort(key=len,reverse=True)
        reps_dict = find_representatives(cliqs, d_l_dict, subg)
        subg = subg.subgraph(reps_dict.keys())
        #merge reps_dict with all_reps_dict
        for key, vals in reps_dict.items():
            if not key in all_reps_dict:
                all_reps_dict[key] = vals
            else:
                #merge represented clusters in a new representative
                newvals = []
                for old_rep in vals:
                    #if statement for bgcs already represented by this 
                    #representative and thus no longer in all_reps_dict
                    if old_rep in all_reps_dict.keys():
                        newv = [v for v in all_reps_dict[old_rep]]
                        newvals += newv
                        del all_reps_dict[old_rep]
                all_reps_dict[key] = set(newvals)
        i+=1
    print("Done. {} representatives chosen for {} bgcs".format(\
        len(all_reps_dict.keys()), g.number_of_nodes()))
    return all_reps_dict

def write_filtered_bgcs(uniq_list, rep_dict, dom_dict, filter_file):
    '''Writes three output files and returns filepath to representatives.csv

    uniq_list: list of strings, bgcs that are not similar to others
    rep_dict: dict of {representative:[represented]}, links representative
        bgcs to bgcs that are filtered out.
    dom_dict: dict of {bgc:[(domains_of_a_gene)]}
    filter_file: str, file path
    Writes three files:
        -filtered_clusterfile.csv: same as clusterfile.csv but without bgcs
        that are filtered out
        -representatives.csv: all the bgcs and their representatives as
        >representative\nbgc1,bgc2\n . also uniq_bgcs are there but just as
        >uniq_bgc1\n>uniq_bgc2\n
        -domstats file only for the representative bgcs
    '''
    rep_file = '{}_representative_bgcs.txt'.format(\
        filter_file.split('_filtered_clusterfile.csv')[0])
    stat_file = '{}_domstats.txt'.format(\
        filter_file.split('_clusterfile.csv')[0])
    domc = Counter()
    with open(filter_file, 'w') as filt, open(rep_file, 'w') as rep:
        for bgc in uniq_list:
            rep.write(">{}\n".format(bgc))
            dom_tups = dom_dict[bgc]
            filt.write("{},{}\n".format(bgc, \
                ','.join(';'.join(gene) for gene in dom_tups)))
            domc.update(dom_tups)
        for bgc in rep_dict.keys():
            rep.write(">{}\n{}\n".format(bgc, ','.join(rep_dict[bgc])))
            dom_tups = dom_dict[bgc]
            filt.write("{},{}\n".format(bgc, \
                ','.join(';'.join(gene) for gene in dom_tups)))
            domc.update(dom_tups)
    with open(stat_file, 'w') as stat:
        stat.write("#Total\t{}\n".format(sum(domc.values())))
        for dom, count in domc.most_common():
            stat.write("{}\t{}\n".format(';'.join(dom),count))
    print("\nFiltered clusterfile containing {} bgcs: {}".format(\
        len(uniq_list)+len(rep_dict.keys()),filt_file))
    print("Representative bgcs file: {}".format(rep_file))
    return rep_file

def remove_infr_doms(clusdict, m_gens, verbose):
    '''Returns clusdict with genes replaced  with (-) if they occur < 3

    clusdict: dict of {cluster:[(domains_of_a_gene)]}
    m_gens: int, minimal distinct genes a cluster must have to be included
    verbose: bool, if True print additional info

    Deletes clusters with 1 unique gene
    '''
    print('\nRemoving domain combinations that occur less than 3 times')
    domcounter = Counter()
    domcounter.update([v for vals in clusdict.values() for v in vals \
        if not v == ('-',)])
    deldoms = [key for key in domcounter if domcounter[key] <= 2]
    print('  {} domain combinations are left, {} are removed'.format(\
        len(domcounter.keys())-len(deldoms),len(deldoms)))
    clus_no_deldoms = {}
    for k,v in clusdict.items():
        newv = [('-',) if dom in deldoms else dom for dom in v]
        doml = len({v for v in newv if not v == ('-',)})
        if doml >= m_gens:
            clus_no_deldoms[k] = newv
        else:
            if verbose:
                print('  {} removed as it has less than min_genes'.format(k))
    print(' {} clusters have less than {} domains and are excluded'.format(\
        len(clusdict.keys()) - len(clus_no_deldoms), m_gens))
    return clus_no_deldoms

def remove_dupl_doms(cluster):
    '''
    Replaces duplicate domains in a cluster with '-', writes domain at the end

    cluster: list of tuples, tuples contain str domain names
    '''
    domc = Counter(cluster)
    dupl = [dom for dom in domc if domc[dom] > 1 if not dom == ('-',)]
    if dupl:
        newclus = [('-',) if dom in dupl else dom for dom in cluster]
        for dom in dupl:
            newclus += [('-',),dom]
    else:
        newclus = cluster
    return newclus

def count_adj(counts, cluster):
    '''Counts all adjacency interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    cluster: list of tuples, genes with domains
    '''
    if len(cluster) == 1:
        return
    for i, dom in enumerate(cluster):
        if i == 0:
            edge = 1
            adj = [cluster[1]]
        elif i == len(cluster)-1:
            edge = 1
            adj = [cluster[i-1]]
        else:
            edge = 2
            adj = [cluster[i-1],cluster[i+1]]
            if adj[0] == adj[1] and adj[0] != ('-',):
                #B2 and N2 counts
                prevdom = cluster[i-1]
                counts[prevdom]['N1'] -= 2
                counts[prevdom]['N2'] += 1
                if dom != ('-',) and dom != prevdom:
                    counts[prevdom]['B1'][dom] -= 2
                    try:
                        counts[prevdom]['B2'][dom] += 1
                    except TypeError:
                        counts[prevdom]['B2'][dom] = 1
        if not dom == ('-',):
            counts[dom]['count'] += 1
            counts[dom]['N1'] += edge
            for ad in adj:
                if ad != ('-',) and ad != dom:
                    try:
                        counts[dom]['B1'][ad] += 1
                    except TypeError:
                        counts[dom]['B1'][ad] = 1

def count_coloc(counts, cluster):
    '''Counts all colocalisation interactions between domains in a cluster

    counts: nested dict { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    cluster: list of tuples, genes with domains
    verbose: bool, if True print additional info
    '''
    N1 = len(cluster)-1
    for dom in cluster:
        if not dom == ('-',):
            counts[dom]['count'] += 1
            counts[dom]['N1'] += N1
            coloc = set(cluster)
            try:
                coloc.remove(('-',))
            except KeyError:
                pass
            coloc.remove(dom)
            for colo in coloc:
                try:
                    counts[dom]['B1'][colo] += 1
                except TypeError:
                    counts[dom]['B1'][colo] = 1

def makehash():
    '''Function to initialise nested dict
    '''
    return defaultdict(makehash)

def count_interactions(clusdict, verbose):
    '''Count all adj and coloc interactions between all domains in clusdict

    clusdict: dict of {cluster:[(gene_with_domains)]}
    verbose: bool, if True print additional info
    Returns two dicts, one dict with adj counts and one with coloc counts
    adj counts:
        { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    coloc counts:
        { dom1:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    '''
    print('\nCounting colocalisation and adjacency interactions')
    all_doms = {v for val in clusdict.values() for v in val}
    all_doms.remove(('-',))

    #initialising count dicts
    adj_counts = makehash()
    for d in all_doms:
        for v in ['count','N1','N2']:
            adj_counts[d][v] = 0
        for w in ['B1','B2']:
            adj_counts[d][w] = makehash()
        #N1: positions adj to one domA, N2: positions adj to two domA
        #B1: amount of domB adj to one domA, B2: positions adj to two domA

    coloc_counts = makehash()
    for d in all_doms:
        for v in ['count','N1']:
            coloc_counts[d][v] = 0
        coloc_counts[d]['B1'] = makehash()
        #N1: all possible coloc positions in a cluster, cluster lenght - 1
        #B1: amount of domB coloc with domA

    for clus in clusdict.values():
        count_adj(adj_counts, clus)
        filt_clus = remove_dupl_doms(clus)
        count_coloc(coloc_counts, filt_clus)
    return(adj_counts, coloc_counts)

def calc_adj_pval_wrapper(count_dict, clusdict, cores, verbose):
    '''Returns list of tuples of corrected pvals for each gene pair

    counts: nested dict { dom1:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    clusdict: dict of {cluster:[(domains_in_a_gene)]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    '''
    print('Calculating adjacency pvalues')
    N = sum([len(values) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=20)
    pvals_ori = pool.map(partial(calc_adj_pval, counts=count_dict, Nall=N), \
        count_dict.items())
    #remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    #to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0],tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0],p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0],p[1])] != 2]
            print('  error with domain pairs {}'.format(', '.join(p_excl)))
            print('  these are excluded')
    #Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method='fdr_by')[1]
    #adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in \
        zip(pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]):
        assert(ab1[0]==ab2[0] and ab1[1]==ab2[1])
        pmax = max([p1,p2])
        ptups.append(((ab1[0],ab1[1]),pmax))
    return ptups

def calc_adj_pval(domval_pair, counts, Nall):
    '''Returns a list of sorted tuples (domA,domB,pval)

    domval_pair: tuple of (domA, {count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w}} )
    Nall: int, all possible positions
    counts: nested dict { domA:{ count:x,N1:y,N2:z,B1:{dom2:v},B2:{dom2:w} } }
    '''
    domA, vals = domval_pair
    #domains without interactions do not end up in pvals
    if not vals['B1'] and not vals['B2']:
        return
    pvals = []
    count = vals['count']
    Ntot = Nall - count
    N1 = vals['N1']
    N2 = vals['N2']
    N0 = Ntot - N1 - N2
    interactions = vals['B1'].keys() | vals['B2'].keys()
    for domB in interactions:
        if domB not in vals['B2']:
            B1 = vals['B1'][domB]
            Btot = counts[domB]['count']
            pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N1, d) \
                for d in range(B1)]) / ncr(Ntot,Btot))
        elif vals['B1'][domB] == 0:
            B2 = vals['B2'][domB]
            Btot = counts[domB]['count']
            pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N2, d) \
                for d in range(B2)]) / ncr(Ntot,Btot))
        else:
            B1 = vals['B1'][domB]
            B2 = vals['B2'][domB]
            Btot = counts[domB]['count']
            pval = float(\
                1 - sum([ncr(N0,Btot-d1-d2) * ncr(N1,d1) * ncr(N2,d2) \
                for d1,d2 in product(range(B1+1),range(B2+1)) \
                if d1+d2 != B1+B2]) / ncr(Ntot,Btot))
        ab_int = sorted((domA,domB))
        pvals.append((ab_int[0],ab_int[1],pval))
    return pvals

def calc_coloc_pval(domval_pair, counts, Nall):
    '''Returns a list of sorted tuples (domA,domB,pval)
    
    domval_pair: tuple of (domA, { count:x,N1:y,B1:{dom2:v,dom3:w } })
    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    Nall: int, all possible positions in all clusters
    '''
    domA, vals = domval_pair
    #domains without interactions do not end up in pvals
    if not vals['B1']:
        return
    pvals = []
    count = vals['count']
    Ntot = Nall - count
    N1 = vals['N1']
    N0 = Ntot - N1
    interactions = vals['B1'].keys()
    for domB in interactions:
        B1 = vals['B1'][domB]
        Btot = counts[domB]['count']
        pval = float(1 - sum([ncr(N0,(Btot-d)) * ncr(N1, d) \
            for d in range(B1)]) / ncr(Ntot,Btot))
        ab_int = sorted((domA,domB))
        pvals.append((ab_int[0],ab_int[1],pval))
    return pvals

def calc_coloc_pval_wrapper(count_dict, clusdict, cores, verbose):
    '''Returns list of tuples of corrected pvals for each domain pair

    counts: nested dict { domA:{ count:x,N1:y,B1:{dom2:v,dom3:w } } }
    clusdict: dict of {cluster:[domains]}
    cores: int, amount of cores to use
    verbose: bool, if True print additional information
    '''
    print('Calculating colocalisation pvalues')
    N = sum([len(remove_dupl_doms(values)) for values in clusdict.values()])
    pool = Pool(cores, maxtasksperchild=50)
    pvals_ori = pool.map(partial(calc_coloc_pval, counts=count_dict, Nall=N), \
        count_dict.items())
    #remove Nones, unlist and sort
    pvals_ori = [lst for lst in pvals_ori if lst]
    pvals_ori = sorted([tup for lst in pvals_ori for tup in lst])
    #to check if there are indeed 2 pvalues for each combination
    check_ps = [(tup[0],tup[1]) for tup in pvals_ori]
    check_c = Counter(check_ps)
    pvals = [p for p in pvals_ori if check_c[(p[0],p[1])] == 2]
    if not len(pvals) == len(pvals_ori):
        if verbose:
            p_excl = [p for p in pvals if check_c[(p[0],p[1])] != 2]
            print('  error with domain pairs {}'.format(', '.join(p_excl)))
            print('  these are excluded')
    #Benjamini-Yekutieli multiple testing correction
    pvals_adj = multipletests(list(zip(*pvals))[2], method='fdr_by')[1]
    #adding adjusted pvals and choosing max
    ptups = []
    for ab1, ab2, p1, p2 in \
        zip(pvals[::2], pvals[1::2], pvals_adj[::2], pvals_adj[1::2]):
        assert(ab1[0]==ab2[0] and ab1[1]==ab2[1])
        pmax = max([p1,p2])
        ptups.append(((ab1[0],ab1[1]),pmax))
    return ptups

def keep_lowest_pval(colocs, adjs):
    '''
    Returns all domain pairs with their lowest pvalue as an edge for nx

    colocs, adjs: list of tuples [((dom1,dom2),pval)]
    Tuples look like (dom1,dom2,{pval:x})
    '''
    pvals = colocs+adjs
    counter = Counter(list(zip(*pvals))[0])
    dupl = sorted([tup for tup in pvals if counter[tup[0]] == 2])
    uniques = [tup for tup in pvals if counter[tup[0]] == 1]
    lowest = []
    for p1,p2 in zip(dupl[::2],dupl[1::2]):
        pmin = min([p1[1],p2[1]])
        lowest.append((p1[0][0],p1[0][1],{'pval':pmin}))
    uniques = [(tup[0][0],tup[0][1],{'pval':tup[1]}) for tup in uniques]
    return lowest+uniques

def visualise_graph(graph, subgraph_list = None, groups = True):
    '''Plots a graph with possible subgraphs in different colours

    graph: networkx graph
    subgraph_list: list of lists of node names that should be coloured
        differently, default = None
    groups: bool, are there groups in the subgraph_list that you want to
        colour differently (True)? or are nodes in subgraph list one seperate
        group (False)
    '''
    cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    options = {'node_size': 2,'width': 0.2}
    pos = nx.spring_layout(graph)
    plt.figure()
    nx.draw_networkx(graph, pos=pos, with_labels = False, node_color='black',\
        **options)
    if subgraph_list:
        if groups:
            for sub in subgraph_list:
                nx.draw_networkx_nodes(graph, pos=pos, nodelist=sub, \
                    node_color=random.choice(cols), **options)
        else:
            nx.draw_networkx_nodes(graph, pos=pos, nodelist=subgraph_list, \
                    node_color='#91bfdb', marker='s', **options)
    plt.show()

def generate_modules_wrapper(pval_edges, sign_cutoff, cores, \
    verbose):
    '''
    Returns a dict with all modules {(module):strictest_pval_cutoff}

    pval_edges: list of tuples, [((dom1),(dom2),pval)]
    sign_cutoff: float, pvalue cutoff
    cores: int, number of cores to use
    verbose: bool, if True print additional information
    '''
    print('\nFinding all modules with a pvalue lower than {}'.format(\
        sign_cutoff))
    sign_pvs = [ptup for ptup in pval_edges if ptup[2]['pval'] <= sign_cutoff]
    print('  {} significant domain pair interactions'.format(len(sign_pvs)))
    pv_values = {pv['pval'] for pv in list(zip(*sign_pvs))[2]}
    #watch out if pv_values gets really big, maybe get 100,000 fixed numbers
    #to loop over if there are more than 100,000 pv_values
    print('  looping through {} pvalue cutoffs'.format(len(pv_values)))
    pool = Pool(cores,maxtasksperchild=100)
    modules = pool.imap(partial(generate_modules, dom_pairs=sign_pvs), \
        pv_values, chunksize=250)
    modules_dict = {}
    for p_mods_pair in modules:
        p = list(p_mods_pair)[0]
        mods = p_mods_pair[p]
        for mod in mods:
            try:
                prev_val = modules_dict[mod]
            except KeyError:
                modules_dict[mod] = p
            else:
                if p < prev_val:
                    modules_dict[mod] = p
    print('{} modules detected'.format(len(modules_dict)))    
    return modules_dict

def generate_modules(sign_cutoff, dom_pairs):
    '''
    Returns modules found with a certain cutoff as {sign_cutoff:{(modules)}}

    sign_cutoff: float, cutoff for detecting modules
    dom_pairs: list of tuples, ('dom1', 'dom2', pvalue)
    Modules are all maximal cliques with length > 2
    '''
    edges = (edge for edge in dom_pairs if edge[2]['pval'] <= sign_cutoff)
    mod_graph = generate_graph(edges, False)
    cliqs = nx.algorithms.clique.find_cliques(mod_graph)
    cliqs = {tuple(sorted(clq)) for clq in cliqs if len(clq) > 2}
    return {sign_cutoff:cliqs}

def write_module_file(outfile,modules,bgc_mod_dict = None):
    '''Write modules with info about the modules to outfile

    outfile: string, path
    modules: dict, {(module_tuple):strictest_detection_cutoff}
    bgc_mod_dict: dict of {bgc: [(modules)]}, default=None
    '''
    print('Writing {} modules to {}'.format(len(modules),outfile))
    with open(outfile, 'w') as out:
        if bgc_mod_dict:
            mod_counts = Counter([mod for modlist in bgc_mod_dict.values() \
                for mod in modlist])
            header = ['Number','N_Occurences','N_Genes','N_Domains',\
                'Strictest_detection_cutoff','Module']
            out.write('{}\n'.format('\t'.join(header)))
            for i,pair in enumerate(sorted(modules.items(), \
                key = itemgetter(1))):
                mod,p = pair
                count = mod_counts[mod]
                domlen = sum(len(m) for m in mod)
                info = [i+1,count,len(mod),domlen,p,\
                    ','.join(';'.join(m) for m in mod)]
                out.write('{}\n'.format('\t'.join(map(str,info))))
        else:
            header = ['Length','Strictest_detection_cutoff',\
                'Module']
            out.write('{}\n'.format('\t'.join(header)))
            for i,pair in enumerate(sorted(modules.items(), key = \
                itemgetter(1))):
                mod,p = pair
                out.write('{}\t{}\t{}\n'.format(len(mod),p,','.join(\
                    ';'.join(m) for m in mod)))

def write_bgcs_and_modules(outfile, clusters, bgc_mod_dict, ranked_mods):
    '''Writes a files containing Bgcname\tDomains\tModules

    outfile: string, file path
    clusters: dict linking bgc name to its domains {bgc:[(genes_domains)]}
    bgc_mod_dict: dict linking bgc to its modules {bgc: [(modules)]}
    ranked_mods: dict of {(module): number}
    Before each domain its index is written seperated with a .
    If there are multiple indeces they are separated with a -
    Modules are seperated by '; '
    outfile:
    clustername\tdom1,dom2;dom3,dom1,dom4\t
        0-2.dom1,3.dom4; 1.dom2;dom3,3.dom4\t
        number_dom1_dom4; number_dom2/3_dom4
    '''
    print('\nWriting file that links bgcs to modules at {}'.format(outfile))
    with open(outfile, 'w') as outf:
        for bgc, doms in clusters.items():
            dom_i = {}
            doms = [';'.join(dom) for dom in doms]
            #find numbers of the doms
            for i,dm in enumerate(doms):
                try:
                    dom_i[dm].append(str(i))
                except KeyError:
                    dom_i[dm] = [str(i)]
            dom_str = {dm:'-'.join(inds)+'.'+dm for dm,inds in dom_i.items()}
            mods = bgc_mod_dict[bgc]
            mods_str = '; '.join([','.join([dom_str[';'.join(d)] for d in \
                mod]) for mod in mods])
            mods_nums = '; '.join(map(str,[ranked_mods[mod] for mod in mods]))
            outf.write('{}\t{}\t{}\t{}\n'.format(\
                bgc,','.join(doms),mods_str,mods_nums))

def link_mods2bgc(bgc, doms, modules):
    '''Returns a tuple of (bgc, [(modules)])

    bgc: string, bgc name
    doms: list of tuples, all domain names from the genes in the bgc
    modules: list of tuples of tuples, each tuple contains the genes of a
        module, each gene contains domains
    '''
    modlist = []
    for mod in modules:
        if set(doms).intersection(mod) == set(mod):
            modlist.append(mod)
    return (bgc,sorted(modlist,key=len))

def link_all_mods2bgcs(bgcs, modules, cores):
    '''Returns a dict of {bgc: [(modules)]}

    bgcs: dict of {bgc: [(domains)])
    modules: list of module tuples of gene tuples
    cores: int, amount of cores to use
    '''
    pool = Pool(cores, maxtasksperchild=100)
    bgcs_mod = pool.starmap(partial(link_mods2bgc, modules=modules), \
        bgcs.items())
    bgc_mod_dict = {pair[0]:pair[1] for pair in bgcs_mod}
    return bgc_mod_dict

def remove_infr_mods(bgc_mod_dict, modules_dict):
    '''Returns updated input where modules that occur < 2 times are removed

    bgc_mod_dict: dict of {bgc: [(modules)]}
    modules: dict  of {mod:[info]}
    '''
    new_bgc_dict = deepcopy(bgc_mod_dict)
    new_mods_dict = deepcopy(modules_dict)
    mod_counts = Counter(modules_dict.keys())
    mod_counts.update([mod for modlist in bgc_mod_dict.values() \
        for mod in modlist])
    #I initialised all mods with 1 so the if statements says < 3 instead of 2
    infr_mods = [mod for mod, count in mod_counts.items() if count < 3]
    for infr_mod in infr_mods:
        del new_mods_dict[infr_mod]
    for bgc,mods in new_bgc_dict.items():
        new_bgc_dict[bgc] = [mod for mod in mods if not mod in infr_mods]
    print('\nRemoving {:.1f}% of modules that occur less than twice'.format(\
        (len(modules_dict)-len(new_mods_dict))/len(modules_dict)*100))
    return new_bgc_dict,new_mods_dict

def read_mods_bgcs(modsfile):
    '''Returns dict, {(mod):strictest_pvalue_cutoff}

    modsfile: string, filename
    '''
    with open(mods_file, 'r') as mods_in:
        mods = {}
        mods_in.readline()
        for line in mods_in:
            line = line.strip().split('\t')
            mod = tuple(tuple(gene.split(';')) for gene in \
                line[-1].split(','))
            mods[mod] = line[-2]
    return mods

def filter_out_domains(clusdict,include_list):
    '''Returns the same dict only keeping domains from the include_list

    clusdict: dict of {bgc:[(domains_in_genes)]}
    include_list: list of str, domains to keep
    '''
    newbgcs = {}
    for bgc, genes in clusdict.items():
        newgenes = []
        for gene in genes:
            if gene == ('-',):
                newgenes.append(gene)
            else:
                ngene = []
                for dom in gene:
                    #check if domain is a subPfam
                    m = re.search(r'_c\d+$',dom)
                    if m:
                        if dom[:m.start()] in include_list:
                            ngene.append(dom)
                    else:
                        if dom in include_list:
                            ngene.append(dom)
                #if gene becomes empty print -
                if not ngene:
                    newgenes.append(('-',))
                else:
                    newgenes.append(tuple(ngene))
        newbgcs[bgc] = newgenes
    return newbgcs

def read_txt(in_file):
    '''Reads text file into list

    in_file: str, file path
    '''
    with open(in_file, 'r') as inf:
        lines = [line.strip() for line in inf]
    return lines

if __name__ == "__main__":
    start = time.time()
    cmd = get_commands()

    #converting genes in each bgc to a combination of domains
    if cmd.start_from_clusterfile:
        if not os.path.isdir(cmd.out_folder):
            f_command = 'mkdir {}'.format(cmd.out_folder)
            subprocess.check_call(f_command,shell=True)
        filepre = os.path.split(cmd.start_from_clusterfile)[-1].split(\
            '.csv')[0]
        clus_file = os.path.join(cmd.out_folder, filepre+'_clusterfile.csv')
        c_command = 'cp {} {}'.format(cmd.start_from_clusterfile,clus_file)
        subprocess.check_call(c_command, shell=True)
    else:
        fasta_folder, exist_fastas = process_gbks(cmd.in_folder, \
            cmd.out_folder, cmd.exclude, cmd.exclude_contig_edge, \
            cmd.min_genes, cmd.cores, cmd.verbose, cmd.use_fastas)
        dom_folder, exist_doms = hmmscan_wrapper(fasta_folder, cmd.hmm_path,\
            cmd.verbose, cmd.cores, exist_fastas, cmd.use_domtabs)
        clus_file = parse_dom_wrapper(dom_folder, cmd.out_folder, \
            cmd.domain_overlap_cutoff, cmd.verbose, exist_doms)

    #filtering clusters based on similarity
    random.seed(595)
    dom_dict = read_clusterfile(clus_file, cmd.min_genes, \
        cmd.verbose)
    doml_dict = {bgc: sum(len(g) for g in genes if not g == ('-',)) \
        for bgc,genes in dom_dict.items()}
    filt_file = '{}_filtered_clusterfile.csv'.format(\
        clus_file.split('_clusterfile.csv')[0])
    similar_bgcs = generate_edges(dom_dict, cmd.sim_cutoff,\
        cmd.cores)
    graph = generate_graph(similar_bgcs, True)
    uniq_bgcs = [clus for clus in dom_dict.keys() if not clus in graph.nodes()]
    all_reps = find_all_representatives(doml_dict, graph)
    if cmd.include_list:
        include_list = read_txt(cmd.include_list)
        dom_dict = filter_out_domains(dom_dict, include_list)
    all_reps_file = write_filtered_bgcs(uniq_bgcs, all_reps, \
        dom_dict, filt_file)

    #detecting modules with statistical approach
    f_clus_dict = read_clusterfile(filt_file, cmd.min_genes, cmd.verbose)
    f_clus_dict_rem = remove_infr_doms(f_clus_dict, cmd.min_genes,cmd.verbose)
    adj_counts, c_counts = count_interactions(f_clus_dict_rem, cmd.verbose)
    adj_pvals = calc_adj_pval_wrapper(adj_counts, f_clus_dict_rem, cmd.cores,\
        cmd.verbose)
    col_pvals = calc_coloc_pval_wrapper(c_counts, f_clus_dict_rem, cmd.cores,\
        cmd.verbose)
    pvals = keep_lowest_pval(col_pvals,adj_pvals)
    mods = generate_modules_wrapper(pvals, cmd.pval_cutoff, cmd.cores,\
        cmd.verbose)
    mod_file = '{}_modules.txt'.format(\
        filt_file.split('_filtered_clusterfile.csv')[0])
    write_module_file(mod_file, mods)
    #linking modules to bgcs and filtering mods that occur less than twice
    bgcs_with_mods_ori = link_all_mods2bgcs(f_clus_dict_rem, mods, cmd.cores)
    bgcs_with_mods, modules = remove_infr_mods(bgcs_with_mods_ori, mods)
    mod_file_f = '{}_filtered_modules.txt'.format(\
        filt_file.split('_filtered_clusterfile.csv')[0])
    write_module_file(mod_file_f,modules,bgcs_with_mods)
    bgcmodfile = '{}_bgcs_with_domains.txt'.format(\
        mod_file.split('_modules.txt')[0])
    rank_mods = {pair[0]:i+1 for i,pair in enumerate(sorted(modules.items(),\
        key=itemgetter(1)))}
    write_bgcs_and_modules(bgcmodfile,f_clus_dict_rem,bgcs_with_mods,\
        rank_mods)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
