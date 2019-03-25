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

'''

from sys import argv
import os
from glob import glob
import subprocess
from Bio import SeqIO
from collections import OrderedDict

def process_gbks(input_folder, exclude, exclude_contig_edge,
    min_genes):
    '''Convert gbk files from input folder to fasta files for each gbk file

    input_folder: str
    exclude: list of str, files will be excluded if part of the file name
        is present in this list
    '''
    if input_folder.endswith('/'):
        out_folder_base, inf = os.path.split(input_folder[:-1])
    else:
        out_folder_base, inf = os.path.split(input_folder)
    out_fasta = os.path.join(out_folder_base, inf+'_fasta')
    if not os.path.isdir(out_fasta):
        subprocess.check_call("mkdir {}".format(out_fasta), shell = True)

    files = glob(os.path.join(input_folder, "*.gbk"))
    processed = 0
    excluded = 0
    for i, file_path in enumerate(files):
        file_name = os.path.split(file_path)[1]
        if any([word in file_name for word in exclude]):
            excluded += 1
            continue
        else:
            convert_gbk2fasta(file_path, out_fasta, exclude_contig_edge,
                min_genes)
            exit()
        processed +=1
    print("Processed {} files and excluded {} files containing {}".format(\
        processed, excluded, ' or '.join(exclude)))

def convert_gbk2fasta(file_path, out_folder, exclude_contig_edge, min_genes):
    '''Convert one gbk file to a fasta file in out_folder

    file_path, out_folder: strings
    '''
    file_name = os.path.split(file_path)[1]
    seqs = OrderedDict()
    try:
        record = next(SeqIO.parse(file_path, 'genbank'))
    except ValueError as e:
        print(" {} in file {}. File will be excluded.".format(e, file_path))
        return
    for feature in record.features:
        if feature.type == 'cluster':
            if "contig_edge" in feature.qualifiers:
                if feature.qualifiers["contig_edge"][0] == "True":
                    if exclude_contig_edge:
                        print(" Contig edge detected in {}".format(file_name))
                        return

        if feature.type == 'CDS':
            print(feature.qualifiers)


if __name__ == "__main__":
    in_folder = argv[1]

    process_gbks(in_folder, exclude = ['final'], exclude_contig_edge = True,
        min_genes = 5)
