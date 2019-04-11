#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of domains and filter the domains based
on similarity.

Usage:
python3 preprocessing_bgcs.py -h

Example usage:
python3 bgc_to_pfam.py -i ../testdata -o ../testdata_domains --hmm_path 
    ../domains/Pfam_100subs_tc.hmm --exclude final -c 12 -e True

Notes:
Only handles gbk files with one cluster

Layout:
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

Required:
python 3.6
hmmscan
Biopython
networkx
'''
import argparse
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter
from functools import partial
from glob import glob, iglob
from itertools import combinations
from multiprocessing import Pool, cpu_count
import networkx as nx
import os
import random
import subprocess

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn bgcs in \
        gbk files into strings of domains using a domain hmm database and to \
        reduce redundancy by filtering out similar bgcs.")
    parser.add_argument("-i", "--in_folder", dest="in_folder", help="Input \
        directory of gbk files", required=True)
    parser.add_argument("--exclude", dest="exclude", default="final",
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
    parser.add_argument("--min_doms", dest="min_doms", default=2,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 2 domains", type=int)
    parser.add_argument("--sim_cutoff", dest="sim_cutoff", default=0.95,
        help="Cutoff for cluster similarity in redundancy filtering (default:\
        0.95)", type=float)
    return parser.parse_args()

def process_gbks(input_folder, output_folder, exclude, exclude_contig_edge,\
    min_genes, cores, verbose):
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
            exclude_contig_edge, min_genes, exclude, verbose), \
            callback=lambda x: done.append(x))
    pool.close()
    pool.join()
    processed = len([val for val in done if val])
    excluded = len([val for val in done if val == False])
    filtered = len([val for val in done if val == None])
    print("Processed {} gbk files into {} fasta files.".format(\
        processed+excluded+filtered, processed))
    print(" excluded {} files containing {}".format(excluded,\
        ' or '.join(exclude)))
    print(" filtered {} files that didn't pass constraints".format(\
        filtered))
    return(out_fasta)

def convert_gbk2fasta(file_path, out_folder, exclude_contig_edge, min_genes,\
    exclude, verbose):
    '''Convert one gbk file to a fasta file in out_folder

    file_path, out_folder: strings
    exclude_contig_edge: bool
    min_genes: int
    verbose: bool, print additional info to stdout

    Returns True for a successful conversion to fasta, None if there is a
    contig edge or min_genes is not passed. False is returned if any of
    the exclude list is in the filename
    '''
    file_name = os.path.split(file_path)[1]
    if any([word in file_name for word in exclude]):
        return False
    name = file_name.strip('.gbk')
    outfile = os.path.join(out_folder, '{}.fasta'.format(name))
    seqs = OrderedDict()
    num_genes = 0
    if not os.path.exists(outfile):
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
                header = ">{}_{}".format(name, num_genes+1)
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
                out.write('{}\n{}\n'.format(seq, seqs[seq]))
    return True

def run_hmmscan(fasta_file, hmm_file, out_folder, verbose):
    """
    Runs hmmscan on fasta file to generate a domtable file

    fasta_file, hmm_file, out_folder: strings of file paths
    verbose: bool
    """
    if os.path.isfile(fasta_file):
        name = os.path.split(fasta_file)[1].split('.fasta')[0]
        out_name = os.path.join(out_folder, name+".domtable")
        log = os.path.join(out_folder, 'hmmlog.txt')
        if not os.path.isfile(out_name):
            hmmscan_cmd = (\
                "hmmscan -o {} --cpu 0 --domtblout {} --cut_tc {} {}".format(\
                log, out_name, hmm_file, fasta_file))
            if verbose:
                print("  " + hmmscan_cmd)
            subprocess.check_call(hmmscan_cmd, shell=True)
        elif verbose:
            print("  {} existed. hmmscan not run again".format(out_name))
    else:
        raise SystemExit("Error running hmmscan: {} doesn't exist".format(\
            fasta_file))

def hmmscan_wrapper(input_folder, hmm_file, verbose, cores):
    '''Runs hmmscan on all fasta files in input_folder hmm_file as hmm db

    fasta_folder, hmm_file: strings of file paths
    verbose: bool
    '''
    if input_folder.endswith('/'):
        out_folder = input_folder[:-7]+'_domtables'
    else:
        out_folder = input_folder[:-6]+'_domtables'
    if not os.path.isdir(out_folder):
        subprocess.check_call("mkdir {}".format(out_folder), shell = True)
    print("\nRunning hmmscan on fastas to generate domtables.")
    files = iglob(os.path.join(input_folder, "*.fasta"))
    pool = Pool(cores, maxtasksperchild=1)
    #maxtasksperchild=1:process respawns after completing 1 task
    for processed, file_path in enumerate(files):
        #run_hmmscan(file_path, hmm_file, out_folder, verbose)
        pool.apply_async(run_hmmscan,args=(file_path, hmm_file,
            out_folder, verbose))
    pool.close()
    pool.join() #make the code in main wait for the pool processes to finish
    print("Processed {} fasta files into domtables.".format(\
        processed+1))
    return out_folder

def parse_domtab(domfile, clus_file, min_overlap, verbose):
    '''Parses domtab into a cluster domain file (csv)

    domfile: string, file path
    clus_file: open file for writing
    min_overlap : float, the amount of overlap two domains must have for it
        to be considered overlap

    clus_file will look like this:
    Clus1,dom1,dom2,-(gene without domain)\\nClus2,dom1..
    '''
    if verbose:
        print("  parsing domtable {}".format(domfile))
    queries = SearchIO.parse(domfile, 'hmmscan3-domtab')
    cds_before = 0
    cluster_doms = [] #domain list for the cluster
    for query in queries:
        dom_matches = []
        cds_num = int(query.id.split('_')[-1])
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
                    if sign_overlap(dom_matches[i][1],dom_matches[j][1],
                        min_overlap):
                        if dom_matches[i][2] >=dom_matches[j][2]:
                            dels.append(j)
                        else:
                            dels.append(i)
        cds_doms = [dom_matches[i][0] for i in range(len(query)) \
            if i not in dels]
        #if a cds has no domains print '-' in output
        gene_gap = cds_num - cds_before -1
        if gene_gap > 0:
            cds_doms = ['-']*gene_gap + cds_doms
        cluster_doms += cds_doms
        cds_before = cds_num
    clus_file.write('{},{}\n'.format(\
        os.path.split(domfile)[-1].split('.domtable')[0],
        ','.join(cluster_doms)))
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

def parse_dom_wrapper(in_folder, out_folder, cutoff, verbose):
    '''Calls parse_domtab on all domtable files to create a clusterfile

    in_folder, out_folder: strings, filepaths
    cutoff: float, cutoff value for domain overlap
    '''
    print("\nParsing domtables from folder {}".format(in_folder))
    domtables = iglob(os.path.join(in_folder, '*.domtable'))
    in_name = os.path.split(in_folder)[1].split('_domtables')[0]
    out_file = os.path.join(out_folder, in_name+'_clusterfile.csv')
    stat_file = os.path.join(out_folder, in_name+'_domstats.txt')
    if not os.path.exists(out_file):
        domc = Counter()
        with open(out_file, 'w') as out:
            for domtable in domtables:
                doms = parse_domtab(domtable, out, cutoff, verbose)
                domc.update(doms)
        with open(stat_file, 'w') as stat:
            stat.write("#Total\t{}".format(sum(domc.values())))
            for dom, count in domc.most_common():
                stat.write("{}\t{}\n".format(dom,count))
    else:
        print("  clusterfile already existed, did not parse again.")
    print("Parsing domtables complete, result in {}".format(out_file))
    print(" statistics about doms in {}".format(stat_file))
    return out_file

def read_clusterfile(infile, m_doms, verbose):
    """Reads a clusterfile into a dictionary

    infile: str, filepath
    m_doms: int, minimum of domains a cluster should have
    verbose: bool, if True print additional info
    clusters with less than m_doms domains are not returned
    """
    print("\nFiltering clusterfile")
    filtered = 0
    with open(infile, 'r') as inf:
        clus_dict = OrderedDict()
        len_dict = OrderedDict()
        for line in inf:
            line = line.strip().split(',')
            clus = line[0]
            doms = line[1:]
            ldoms = len([dom for dom in doms if not dom == '-'])
            if ldoms < m_doms:
                filtered +=1
                if verbose:
                    print("  excluding {} less than min domains".format(clus))
                continue
            if not clus in clus_dict.keys():
                clus_dict[clus] = doms
                len_dict[clus] = ldoms
            else:
                print("Clusternames not unique, {} read twice".format(clus))
    print("Done. Keeping {} clusters".format(len(clus_dict)))
    print(" {} clusters have less than {} domains".format(filtered,m_doms))
    return clus_dict, len_dict

def calc_adj_index(clus1, clus2):
    '''Returns the adjacency index between two clusters

    clus1, clus2: list of strings, domainlist of a cluster

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

    clus1, clus2: list of strings, domainlist of a cluster
    '''
    one_in_two = all([dom in clus2 for dom in clus1 if not dom == '-'])
    two_in_one = all([dom in clus1 for dom in clus2 if not dom == '-'])
    if one_in_two or two_in_one:
        return True
    return False

def generate_edges(nodes, dom_dict, cutoff, cores):
    '''Returns a pair of clusters in a tuple if ai/contained above cutoff

    nodes: list of strings, clusternames
    dom_dict: dict {clus1:[domains]}, clusters linked to domains
    cutoff: float, between 0-1, when clusters are similar
    cores: int, amount of cores used for calculation
    '''
    print("\nGenerating similarity scores")
    pairs = combinations(clus_names, 2)
    pool = Pool(cores, maxtasksperchild = 100)
    #I could add imap if this is still too slow for antismashdb
    #with imap specify chunksize as number of pairs?
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
    contained = is_contained(d_dict[p1], d_dict[p2])
    ai = calc_adj_index(d_dict[p1],d_dict[p2])
    if ai == None:
        print(ai,pair)
    if contained or ai > cutoff:
        # print(pair,ai,contained)
        return(p1,p2,{'ai':ai,'contained':contained})

def generate_graph(edges):
    '''Returns a networkx graph

    edges: list of tuples, (pair1,pair2,{attributes})
    '''
    g = nx.Graph()
    g.add_edges_from(edges)
    print('Generated graph with:')
    print(' {} nodes'.format(g.number_of_nodes()))
    print(' {} edges'.format(g.number_of_edges()))
    return g

def find_representatives(clqs, d_l_dict, graph):
    '''
    Returns {representative:[clique]} based on bgc with most domains in clique

    clqs: list of lists of strings, cliques of clusters
    d_l_dict: dict of {clus_name:amount_of_domains(int)}
    graph: networkx graph structure of the cliques
    Returns also the representative in the clique list (dict values)
    The longest cluster is chosen (most domains). If there are multiple
        longest clusters then the cluster with the least connections is
        chosen (to preserve most information).
    '''
    reps_dict = OrderedDict()
    dels = set() #set of nodes for which a representative has been found
    clqs = sorted(clqs, key=len, reverse=True)
    for cliq in clqs:
        cliq = [clus for clus in cliq if not clus in dels]
        if cliq:
            domlist = [(clus,d_l_dict[clus]) for clus in cliq]
            maxdoml = max([doms[1] for doms in domlist])
            clus_maxlen = [clus for clus, doml in domlist \
                if doml == maxdoml]
            if len(clus_maxlen) > 1:
                min_degr = min([deg for clus, deg in graph.degree(clus_maxlen)])
                random.seed(1)
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
    graph: networkx graph structure containing the cliques
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
        cliqs = sorted(cliqs, key=len)
        cliqs = [cl for cl in cliqs if len(cl) > 1]
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
    '''Returns filepaths to filtered_clusterfile.csv and representatives.csv

    uniq_list: list of strings, bgcs that are not similar to others
    rep_dict: dict of {representative:[represented]}, links representative
        bgcs to bgcs that are filtered out.
    Writes two files:
        -filtered_clusterfile.csv: same as clusterfile.csv but without bgcs
        that are filtered out
        -representatives.csv: all the bgcs and their representatives as
        >representative\nbgc1,bgc2\n . also uniq_bgcs are there but just as
        >uniq_bgc1\n>uniq_bgc2\n
    '''
    rep_file = '{}_representative_bgcs.txt'.format(\
        filter_file.split('_filtered_clusterfile.csv')[0])
    with open(filter_file, 'w') as filt, open(rep_file, 'w') as rep:
        for bgc in uniq_list:
            rep.write(">{}\n".format(bgc))
            filt.write("{},{}\n".format(bgc, ','.join(dom_dict[bgc])))
        for bgc in rep_dict.keys():
            rep.write(">{}\n{}\n".format(bgc, ','.join(rep_dict[bgc])))
            filt.write("{},{}\n".format(bgc, ','.join(dom_dict[bgc])))
    print("\nFiltered clusterfile containing {} bgcs: {}".format(\
        len(uniq_list)+len(rep_dict.keys()),filt_file))
    print("Representative bgcs file: {}".format(rep_file))
    return rep_file

if __name__ == "__main__":
    cmd = get_commands()

    fasta_folder = process_gbks(cmd.in_folder, cmd.out_folder, cmd.exclude,
        cmd.exclude_contig_edge, cmd.min_genes, cmd.cores, cmd.verbose)
    dom_folder = hmmscan_wrapper(fasta_folder, cmd.hmm_path, cmd.verbose,
        cmd.cores)
    clus_file = parse_dom_wrapper(dom_folder, cmd.out_folder, \
        cmd.domain_overlap_cutoff, cmd.verbose)

    random.seed(1)
    dom_dict, doml_dict = read_clusterfile(clus_file, cmd.min_doms, \
        cmd.verbose)
    filt_file = '{}_filtered_clusterfile.csv'.format(\
        clus_file.split('_clusterfile.csv')[0])
    clus_names = list(dom_dict.keys())
    similar_bgcs = generate_edges(clus_names, dom_dict, cmd.sim_cutoff,\
        cmd.cores)
    graph = generate_graph(similar_bgcs)
    uniq_bgcs = [clus for clus in clus_names if not clus in graph.nodes()]
    all_reps = find_all_representatives(doml_dict, graph)
    all_reps_file = write_filtered_bgcs(uniq_bgcs, all_reps, dom_dict, \
        filt_file)
