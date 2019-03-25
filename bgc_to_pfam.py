#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert BGCs into strings of PFAM domains.

Usage:
    python3 bgc_to_pfam.py <input_folder> <pfam_path> <output_folder>

Fasta files created in input_folder_fasta

TODO:
    Keep contig edges into account?
'''

from sys import argv
import os
from glob import glob
import subprocess

def process_gbks(input_folder, exclude = ['final']):
    '''Convert gbk files to fasta files for each gbk files

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

    if os.path.isfile(input_folder):
        files = [inputfolder]
    else:
        files = glob(os.path.join(input_folder, "*.gbk"))
    processed = 0
    excluded = 0
    for i, file_path in enumerate(files):
        file_name = os.path.split(file_path)[1]
        if any([word in file_name for word in exclude]):
            excluded += 1
            continue
        else:
            #convert_gbk2fasta(file_name, out_fasta)
            pass
        processed +=1
    print("Processed {} files and exlcuded {} files containing {}".format(\
        processed, excluded, ' or '.join(exclude)))




if __name__ == "__main__":
    in_folder = argv[1]

    process_gbks(in_folder)
