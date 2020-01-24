#!usr/bin/env python3
'''
Author: Joris Louwen

Part of iPRESTO, Bioinformatics group Wageningen University.
PI: Marnix Medema

Script to link NPs known to be produced by certain organism to BGCs through
sub-cluster analysis with iPRESTO.

Required:
-output file of NPatlas for a substructure search
    (https://www.npatlas.org/joomla/index.php/search/advanced-search)
-iPRESTO output (PRESTO-TOP output), sub-cluster motifs. Sub-cluster motifs
    should be grouped by motif. This is for example the 
    matches_per_topic_filtered.txt

Usage:
python3 link_np_to_bgc.py -h
'''

import argparse
from glob import glob, iglob
import os
import time

def get_commands():
    parser = argparse.ArgumentParser(description="Script to link NPs known \
    to be produced by certain organism to BGCs through sub-cluster analysis \
    with iPRESTO. Needed are NPatlas output (tsv) and PRESTO-TOP output.")
    parser.add_argument("-i", "--in_file", help="Input substructure search \
        file downloaded from NPatlas output", required=True)
    parser.add_argument("-o", "--out_file", dest="out_file", 
        required=True, help="Output file")
    parser.add_argument("-m", "--motifs_file", help="File containing\
        sub-cluster motifs, such as matches_per_topic_filtered.txt")
    parser.add_argument("-s", "--subcluster_selection", help="The sub-cluster\
        motif (number) that you are interested in that is annotated with the \
        substructure from the NPatlas search")
    parser.add_argument("-t", "--taxonomy_file", help="File linking BGCs\
        to organism names, tsv of bgc\torganism")
    return parser.parse_args()

def extract_subcl_motif(infile, motif):
    '''Extract all BGCs with motif_num as {bgc: [info_motif_match]}

    infile: str, filepath to motif file
    motif: str, name of motif
    '''
    right_motif = False
    bgc_dict = {}
    print(infile)
    with open(infile, 'r') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('#'):
                right_motif = False
                #header of motif matches
                motif_in_file = line.split(', ')[0].split(' ')[-1]
                if motif_in_file == motif:
                    right_motif = True
            else:
                if right_motif:
                    line = line.split('\t')
                    bgc = line.pop(3)
                    bgc_dict[bgc] = line #add filter??
    return bgc_dict


if __name__ == "__main__":
    start = time.time()
    print("\nStart")
    cmd = get_commands()

    bgcs_in_motif = extract_subcl_motif(cmd.motifs_file,
        cmd.subcluster_selection)

    print(bgcs_in_motif, len(bgcs_in_motif))

    end = time.time()
    t = end-start
    t_str = '{}h{}m{}s'.format(int(t/3600),int(t%3600/60),int(t%3600%60))
    print('\nScript completed in {}'.format(t_str))
