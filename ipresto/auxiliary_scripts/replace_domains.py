#!/usr/bin/env python3
'''
Author: Joris Louwen
Student number: 960516530090

Script to delete certain hmm models from a hmm flatfile based on
domain names provided in a text file by creating a new hmm flatfile
without said domains.

Usage:
    python3 replace_domains.py <main_hmm> <delete_domains.txt> <out_hmm>

Notes:
domains in delete_domains should be on seperate lines and match exactly
'''

from sys import argv
from itertools import groupby

def read_dels(in_file):
    '''Reads text file into list

    in_file: str, file path
    '''
    with open(in_file, 'r') as inf:
        lines = [line.strip() for line in inf]
    return lines

def delete_entries(in_file, out_file, dels):
    '''Writes a new file with elements from dels deleted in in_file

    in_file, out_file: strings of file paths
    dels: list of str, names to delete from in_file
    '''
    with open(in_file, 'r') as inf, open(out_file, 'w') as outf:
        #returns a bool in boolval and an iterable in lines
        for boolval, lines in groupby(inf, lambda line:\
            line.startswith('HMMER3/f')):
            if boolval:
                header = next(lines) #will only contain 1 line
            else:
                record = [line for line in lines]
                name = record[0].rstrip('\n').split()[1]
                if name in dels:
                    print("Excluding {}".format(name))
                    continue
                else:
                    outf.write(header)
                    outf.write('{}'.format(''.join(record)))

if __name__ == "__main__":
    hmm_file = argv[1]
    del_file = argv[2]
    output = argv[3]

    del_list = read_dels(del_file)
    delete_entries(hmm_file, output, del_list)
