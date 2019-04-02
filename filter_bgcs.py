#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to convert domainlist of bgcs into graphs based on a Adjacency index
and filter out redundant bgcs.

Usage:
python3 bgc_to_pfam.py -h

Example usage:
python3 filter_bgcs.py

Notes:

Layout:

Required:

'''
import os
from glob import glob, iglob
import subprocess
from collections import OrderedDict, Counter
import argparse

def get_commands():
    parser = argparse.ArgumentParser(description="A script to turn a \
        clusterfile with domains into graphs and filter out redundant bgcs")
    parser.add_argument("-i", "--in_file", dest="in_file", help="Input \
        file", required=True)
    parser.add_argument("-o", "--out_folder", dest="out_folder", 
        required=True, help="Output directory, this will contain all output \
        data files.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(), 
        help="Set the number of cores the script may use (default: use all \
        available cores)", type=int)
    parser.add_argument("-v", "--verbose", dest="verbose", required=False,
        action="store_true", default=False, help="Prints more detailed \
        information.")
    parser.add_argument("--min_doms", dest="min_doms", default=5,
        help="The minimum amount of domains in a BGC to be included in the \
        analysis. Default is 0 domains", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    cmd = get_commands()
