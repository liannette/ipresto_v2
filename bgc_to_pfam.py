#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of PFAM domains.

Usage:
    python3 bgc_to_pfam.py <input_folder> <pfam_path> <output_folder>

Notes:
Fasta files created in input_folder_fasta
Only handles gbk files with one cluster
domtables in input_folder_domtables

TODO:
write excluded files to a file and the reason of exclusion (now stdout if verbose).

Layout:
process_gbks
convert_gbk2fasta
run_hmmscan
hmmscan_wrapper
'''

import sys
import os
from glob import glob
import subprocess
from Bio import SeqIO
from collections import OrderedDict
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
    parser.add_argument("-v", "--verbose", dest="verbose", 
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
    out_fasta = os.path.join(output_folder, inf+'_fasta'
    if not os.path.isdir(out_fasta):
        subprocess.check_call("mkdir {}".format(out_fasta), shell = True)
    print("Processing GBK files.")
    files = glob(os.path.join(input_folder, "*.gbk"))
    processed = 0
    excluded = 0
    filtered = 0
    for i, file_path in enumerate(files):
        if processed % 1000 == 0:
            print(" processed {} files".format(processed))
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
        name = os.path.split(fasta_file)[1].strip('.fasta')
        out_name = os.path.join(out_folder, name+".domtable")
        log = os.path.join(out_folder, 'hmmlog.txt')
        hmmscan_cmd = (\
            "hmmscan -o {} --cpu 0 --domtblout {} --cut_tc {} {}".format(\
            log, out_name, hmm_file, fasta_file))
        if verbose:
            print("  " + hmmscan_cmd)
        subprocess.check_call(hmmscan_cmd, shell=True)
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
    print("Running hmmscan on fastas to generate domtables.")
    files = glob(os.path.join(input_folder, "*.fasta"))
    processed = 0
    pool = Pool(cores, maxtasksperchild=1)
    for i, file_path in enumerate(files):
        file_name = os.path.split(file_path)[1]
        #run_hmmscan(file_path, hmm_file, out_folder, verbose)
        pool.apply_async(run_hmmscan,args=(file_path, hmm_file, out_folder, verbose))
        processed +=1
    pool.close()
    pool.join()
    print("Processed {} fasta files into domtables.".format(\
        processed))

if __name__ == "__main__":
    #in_folder = sys.argv[1]
    #hmm = sys.argv[2]
    cmd = get_commands()

    if cmd.in_folder.endswith('/'):
        fasta_folder = in_folder[:-1]+'_fasta'
    else:
        fasta.folder = in_folder+'_fasta'
    #fasta_folder = process_gbks(cmd.in_folder, cmd.exclude,
    #    cmd.exclude_contig_edge, cmd.min_genes, cmd.verbose)
    hmmscan_wrapper(fasta_folder, cmd.hmm, cmd.verbose)
