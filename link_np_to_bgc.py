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
    with iPRESTO.")
    parser.add_argument("-i", "--in_file", help="Input substructure search \
        file downloaded from NPatlas output", required=True)
    parser.add_argument("-o", "--out_file", dest="out_file", 
        required=True, help="Output file")
    parser.add_argument("-m", "--subcluster_motifs", help="File containing\
        sub-cluster motifs, such as matches_per_topic_filtered.txt")
    parser.add_argument("-s", "--subcluster_selection", help="The sub-cluster\
        motif (number) that you are interested in that is annotated with the \
        substructure from the NPatlas search")
    return parser.parse_args()
    

if __name__ == "__main__":
    
