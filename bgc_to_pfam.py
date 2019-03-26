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
    exclude_contig_edge: bool
    min_genes: int
    '''
    if input_folder.endswith('/'):
        out_folder_base, inf = os.path.split(input_folder[:-1])
    else:
        out_folder_base, inf = os.path.split(input_folder)
    out_fasta = os.path.join(out_folder_base, inf+'_fasta')
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
                min_genes)
            if not done:
                filtered +=1
            #if i > 5:
                #exit()
        processed +=1
    print("Processed {} files and excluded {} files containing {}".format(\
        processed, excluded, ' or '.join(exclude)))
    print(" Filtered {} files that didn't pass constraints".format(filtered))
    return(out_fasta)

def convert_gbk2fasta(file_path, out_folder, exclude_contig_edge, min_genes):
    '''Convert one gbk file to a fasta file in out_folder

    file_path, out_folder: strings
    exclude_contig_edge: bool
    min_genes: int
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
                            print(" Excluding {}: {}".format(file_name,\
                                "contig edge"))
                            return
            if feature.type == 'CDS':
                #print(feature.qualifiers)
                header = ">{}_{}".format(name, num_genes+1)
                seqs[header] = feature.qualifiers['translation'][0]
                num_genes +=1

        if num_genes < min_genes:
            print(" Excluding {}: less than {} genes".format(file_path,\
                min_genes))
            return
        with open(outfile, 'w') as out:
            for seq in seqs:
                out.write('{}\n{}\n'.format(seq, seqs[seq]))
    return True

def runHmmScan(fasta_file, hmmPath, outputdir, verbose):
    """ Runs hmmscan command on a fasta file with a single core to generate a
    domtable file"""
    hmmFile = os.path.join(hmmPath,"Pfam-A.hmm")
    if os.path.isfile(fasta_file):
        name = ".".join(fasta_file.split(os.sep)[-1].split(".")[:-1])
        outputName = os.path.join(outputdir, name+".domtable")
        
        hmmscan_cmd = "hmmscan --cpu 0 --domtblout {} --cut_tc {} {}".format(\
            outputName, hmmFile, fasta_file)
        if verbose == True:
            print("   " + hmmscan_cmd)
        subprocess.check_output(hmmscan_cmd, shell=True)
    else:
        raise SystemExit("Error running hmmscan: {} doesn't exist".format(\
            fasta_file))

if __name__ == "__main__":
    in_folder = argv[1]

    fasta_folder = process_gbks(in_folder, exclude = ['final'],
        exclude_contig_edge = True, min_genes = 5)
