#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of domains.

Usage:
python3 bgc_to_pfam.py -h

Example usage:
python3 bgc_to_pfam.py -i ../testdata -o ../testdata_domains --hmm_path 
    ../domains/Pfam_100subs_tc.hmm --exclude final -c 12 -e True

Notes:
Only handles gbk files with one cluster
Multiple cores are only used for hmmscan step

Layout:
process_gbks
convert_gbk2fasta
run_hmmscan
hmmscan_wrapper
parse_domtab
sign_overlap
parse_dom_wrapper

Required:
hmmscan
'''

import sys
import os
from glob import glob, iglob
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from collections import OrderedDict, Counter
from multiprocessing import Pool, cpu_count
import argparse

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn gbk files into \
        strings of domains using a domain hmm database")
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
    return parser.parse_args()

def process_gbks(input_folder, output_folder, exclude, exclude_contig_edge,\
    min_genes, verbose):
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
    processed = 0
    excluded = 0
    filtered = 0
    for file_path in files:
        if processed % 1000 == 0:
            print(" converted {} files".format(processed))
        file_name = os.path.split(file_path)[1]
        if any([word in file_name for word in exclude]):
            excluded += 1
            continue
        else:
            done = convert_gbk2fasta(file_path, out_fasta, exclude_contig_edge,
                min_genes, verbose)
            if not done:
                filtered +=1
        processed +=1
    print("Processed {} gbk files into {} fasta files.".format(\
        processed+excluded, processed-filtered))
    print(" excluded {} files containing {}".format(excluded,\
        ' or '.join(exclude)))
    print(" filtered {} files that didn't pass constraints".format(\
        filtered))
    return(out_fasta)

def convert_gbk2fasta(file_path, out_folder, exclude_contig_edge, min_genes,\
    verbose):
    '''Convert one gbk file to a fasta file in out_folder

    file_path, out_folder: strings
    exclude_contig_edge: bool
    min_genes: int
    verbose: bool, print additional info to stdout
    '''
    file_name = os.path.split(file_path)[1]
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

if __name__ == "__main__":
    cmd = get_commands()

    fasta_folder = process_gbks(cmd.in_folder, cmd.out_folder, cmd.exclude,
        cmd.exclude_contig_edge, cmd.min_genes, cmd.verbose)
    dom_folder = hmmscan_wrapper(fasta_folder, cmd.hmm_path, cmd.verbose,
        cmd.cores)
    parse_dom_wrapper(dom_folder, cmd.out_folder, cmd.domain_overlap_cutoff,
        cmd.verbose)
