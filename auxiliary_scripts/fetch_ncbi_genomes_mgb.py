#!/usr/bin/env python3
'''
Fetch a bunch of genomes specified in a list of plain text blast result files.
Author: Joris Louwen
'''

from Bio import Entrez
from collections import defaultdict
import os
import argparse
import subprocess

def get_commands():
    parser = argparse.ArgumentParser(description="A script to fetch genomes\
        from the ncbi api by looking through blast plain txt result files for\
        hits in organisms above 50% identity.")
    parser.add_argument('-i','--in_file',help='Input files of blast results\
        where organsims are stated within [], and per. ident is the 5th\
        result column seperated by whitespaces. Can be more than one file.',
        nargs='+')
    parser.add_argument('-o','--out_file',help='Directory of output files.')
    return parser.parse_args()

def retrieve_genomes(in_files,out_folder):
    '''
    Fetch all nucleotide entries for organisms in all files, write2out_folder

    in_files: list of str, file locations
    out_folder: str, file location
    '''
    organisms = parse_organisms(in_files)

def parse_organisms(in_files):
    '''Returns list of organism strings if %ID above 50%

    in_files: list of input files
    '''
    orgn_list = []
    for in_file in in_files:
        with open(in_file,'r') as inf:
            start = False
            for line in inf:
                if line.startswith('Description'):
                    start = True
                if start:
                    orgn_raw = line.split('[')[-1]
                    orgn = orgn_raw.split(']')[0]
                    perc_id = float(line.strip().split('\s')[-2])
                    if perc_id > 50:
                        orgn_list.append(orgn)
    print(orgn_list)

if __name__ == '__main__':
    cmd = get_commands()

    if not os.path.isdir(cmd.out_file):
        subprocess.check_call('mkdir {}'.format(cmd.out_file),shell=True)

    retrieve_genomes(cmd.in_file, cmd.out_file)

    '''
    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    >>> handle = Entrez.esearch(db="pubmed", term="biopython")
    >>> record = Entrez.read(handle)
    >>> "19304878" in record["IdList"]
    True
    >>> print(record["IdList"])
    ['28011774', '24929426', '24497503', '24267035', '24194598', ..., '14871861']
    
    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
    >>> handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
    >>> print(handle.read())
    '''
