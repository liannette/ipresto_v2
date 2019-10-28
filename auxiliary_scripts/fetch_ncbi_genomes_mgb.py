#!/usr/bin/env python3
'''
Fetch a bunch of genomes specified in a list of plain text blast result files.
Author: Joris Louwen
'''

from Bio import Entrez
from collections import defaultdict
from glob import glob
import os
import argparse
import subprocess
import re
import time
import urllib

def get_commands():
    parser = argparse.ArgumentParser(description="A script to fetch genomes\
        from the ncbi api by looking through blast plain txt result files for\
        hits in organisms above 50% identity.")
    parser.add_argument('-i','--in_file',help='Input files of blast results\
        where organsims are stated within [], and per. ident is the 5th\
        result column seperated by whitespaces. Can be more than one file.\
        Can also be one directory with files (only one directory)',
        nargs='+')
    parser.add_argument('-s','--strain_list',help='Txt file with a strain\
        name on each line, optional. Overwrites -i',default=False)
    parser.add_argument('-o','--out_file',help='Directory of output files.')
    return parser.parse_args()

def retrieve_genomes(in_files,out_folder,orgn_list=False):
    '''
    Fetch all nucleotide entries for organisms in all files, write2out_folder

    in_files: list of str, file locations
    out_folder: str, file location
    '''
    if orgn_list:
        organisms = orgn_list
    else:
        organisms = parse_organisms(in_files)
    results = []
    try_again = []
    for orgn in organisms:
        try:
            results.append(retrieve_genome(orgn,out_folder))
        except (urllib.error.URLError,RuntimeError):
            try_again.append(orgn)
    for orgn in try_again:
        try:
            results.append(retrieve_genome(orgn,out_folder))
        except (urllib.error.URLError,RuntimeError):
            print('No assembly found for',orgn)
            results.append(None)
    print('{} assemblies downloaded successfully, {} failed'.format(\
        len([1 for a in results if a==True]),len([1 for a in results \
        if a==None])))

def retrieve_genome(orgn, out_folder):
    '''Download assembly of orgn into out_folder

    If no assembly can be found return None
    '''
    print('Downloading assembly for',orgn)
    query = '"{}" AND ("complete genome" OR chromosome OR project)'.format(\
        orgn)
    Entrez.email = 'joris.louwen@wur.nl'
    handle = Entrez.esearch(db="nucleotide", term=query,idtype="acc")
    try:
        record = Entrez.read(handle)
    except IndexError:
        print('No assembly found for',orgn)
        return
    match = None
    for gbk_acc in record['IdList']:
        # gbk_acc = record['IdList'][0]
        handle = Entrez.efetch(db="nucleotide", id=gbk_acc,\
            rettype="gb", retmode="text")
        gbk_data = handle.read()[:1000]
        match = re.search(r'Assembly: (.+?\.\d)',gbk_data)
        if match:
            break
    if not match:
        print('No assembly found for',orgn)
        print(record)
        return
    asmbl = match.group(1)
    ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/'.format(\
        asmbl[:3],asmbl[4:7],asmbl[7:10],asmbl[10:13])
    download_command = (
        'wget -q --recursive -e robots=off --reject "index.html"' +
        ' --no-host-directories --cut-dirs=6 {}'.format(ftp_url) +
        ' -P {}'.format(out_folder))
    out = glob(out_folder+asmbl+'*.gbff')
    if not out:
        subprocess.check_call(download_command,shell=True)
        out_dir = glob(out_folder+asmbl+'*/')
        gzip_genome = glob(os.path.join(out_dir[0],\
            asmbl+'*genomic.gbff.gz'))[0]
        out_genome = os.path.join(out_folder,\
            os.path.split(gzip_genome)[1].split('.gz')[0])
        subprocess.check_call('gunzip -c {} > {}'.format(gzip_genome,out_genome),\
            shell=True)
        subprocess.check_call('rm -r {}'.format(out_dir[0]),shell=True)
    else:
        print('{} existed. Did not download {} again'.format(asmbl,orgn))
    return True

def parse_organisms(in_files):
    '''Returns list of organism strings if %ID above 50%

    in_files: list of input files, xml blast output
    '''
    orgn_set = set()
    for in_file in in_files:
        with open(in_file,'r') as inf:
            rec = False
            for line in inf:
                if not line.strip():
                    continue
                line = line.strip()
                if line.startswith('<Hit_def>'):
                    rec = True
                    matches = re.findall(r'\[[A-Z]\w+\s.+?\]',line)
                    orgns = [m.strip('[]') for m in matches]
                if rec:
                    if line.startswith('<Hsp_identity>'):
                        ident = int(line.split('>')[1].split('<')[0])
                    if line.startswith('<Hsp_align-len>'):
                        al_len = int(line.split('>')[1].split('<')[0])
                        if ident/al_len >= 0.5:
                            rec = False
                            orgn_set.update(orgns)
    return list(orgn_set)

if __name__ == '__main__':
    start = time.time()
    cmd = get_commands()

    if not os.path.isdir(cmd.out_file):
        subprocess.check_call('mkdir {}'.format(cmd.out_file),shell=True)

    if os.path.isdir(cmd.in_file[0]):
        in_files_main = glob(os.path.join(cmd.in_file[0],'*.xml'))
    else:
        in_files_main = cmd.in_file

    retrieve_genomes(in_files_main, cmd.out_file, cmd.strain_list)

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
