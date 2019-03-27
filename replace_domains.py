#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to delete certain hmm models from a hmm flatfile based on
domain names provided in a text file

Usage:
    python3 replace_domains.py <main_hmm> <delete_domains.txt>

Notes:
domains in delete_domains should be on seperate lines
'''

from sys import argv

def read_dels(filename):
    '''Reads text file into list
    '''
    with open(filename, 'r') as infile:
        lines = [line.strip() for line in infile]
    return lines

if __name__ == "__main__":
    hmm_file = argv[1]
    del_file = argv[2]

    del_list = read_dels(del_file)
    print(del_list)
