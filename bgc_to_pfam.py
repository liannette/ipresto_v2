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

TODO:
    write excluded files to a file and the reason of exclusion.

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

def process_gbks(input_folder, exclude, exclude_contig_edge,\
    min_genes, verbose):
    '''Convert gbk files from input folder to fasta files for each gbk file

    input_folder: str
    exclude: list of str, files will be excluded if part of the file name
        is present in this list
    exclude_contig_edge: bool
    min_genes: int
    verbose: bool, print additional info to stdout
    '''
    if input_folder.endswith('/'):
        out_folder_base, inf = os.path.split(input_folder[:-1])
    else:
        out_folder_base, inf = os.path.split(input_folder)
    out_fasta = os.path.join(out_folder_base, inf+'_fasta')
    if not os.path.isdir(out_fasta):
        subprocess.check_call("mkdir {}".format(out_fasta), shell = True)
    sys.stderr.write("Processing GBK files.")
    files = glob(os.path.join(input_folder, "*.gbk"))
    processed = 0
    excluded = 0
    filtered = 0
    for i, file_path in enumerate(files):
        if processed % 1000 == 0:
            sys.stderr.write(" processed {} files".format(processed))
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
    sys.stderr.write("Processed {} gbk files into {} fasta files.".format(\
        processed+excluded, processed-filtered))
    sys.stderr.write(" excluded {} files containing {}".format(excluded,\
        ' or '.join(exclude)))
    sys.stderr.write(" filtered {} files that didn't pass constraints".format(\
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
    Runs hmmscan on fasta file with a single core to generate a domtable file

    fasta_file, hmm_file, out_folder: strings of file paths
    verbose: bool
    """
    if os.path.isfile(fasta_file):
        name = os.path.split(file_path)[1].strip('.fasta')
        out_name = os.path.join(out_folder, name+".domtable")
        hmmscan_cmd = "hmmscan --cpu 0 --domtblout {} --cut_tc {} {}".format(\
            out_name, hmm_file, fasta_file)
        if verbose:
            print("  " + hmmscan_cmd)
        subprocess.check_output(hmmscan_cmd, shell=True)
    else:
        raise SystemExit("Error running hmmscan: {} doesn't exist".format(\
            fasta_file))

def hmmscan_wrapper(input_folder, hmm_file, verbose):
    '''Runs hmmscan on all fasta files in input_folder hmm_file as hmm db

    fasta_folder, hmm_file: strings of file paths
    '''
    if input_folder.endswith('/'):
        out_folder_base, inf = os.path.split(input_folder[:-1])
    else:
        out_folder_base, inf = os.path.split(input_folder)
    out_folder = os.path.join(out_folder_base, inf+'_domtables')
    if not os.path.isdir(out_folder):
        subprocess.check_call("mkdir {}".format(out_folder), shell = True)
    sys.stderr.write("Running hmmscan on fastas to generate domtables.")
    files = glob(os.path.join(input_folder, "*.fasta"))
    processed = 0
    for i, file_path in enumerate(files):
        if processed % 1000 == 0:
            sys.stderr.write(" processed {} files".format(processed))
        file_name = os.path.split(file_path)[1]
        run_hmmscan(file_path, hmm_file, out_folder, verbose)
        processed +=1
    sys.stderr.write("Processed {} fasta files into domtables.".format(\
        processed))

if __name__ == "__main__":
    in_folder = sys.argv[1]

    verbose = False
    fasta_folder = process_gbks(in_folder, exclude = ['final'],
        exclude_contig_edge = True, min_genes = 5, verbose)
