#!/usr/bin/env python3
'''
Author: Joris
Quickly check stats on sub-cluster motifs with BGCs, using a filtered list of
sub-cluster motifs to use.
'''

from collections import Counter
import matplotlib.pyplot as plt
from sys import argv



def read_txt_to_list(filename):
    '''
    Read txt file with fasta like headers to {header:{class:'',motif:[info]}}

    filename: str, filepath of input file
    '''
    bgc_subcl = {}
    with open(filename, 'r') as inf:
        for line in inf:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                header = line[1:]
                bgc_subcl[header] = {}
            else:
                if line.startswith('class'):
                    bgc_subcl[header]['class'] = line.split('=')[-1]
                else:
                    line = line.split('\t')
                    bgc_subcl[header][line[0]] = line[1:]
    return bgc_subcl

def read_filter_subcl_mots(filename):
    '''Read txt tab delim txt file, return only first column to list of str

    filename: str, filepath of input file
    '''
    subcl_filt = []
    header_tokens = '#ST'
    with open(filename, 'r', encoding="utf-8") as inf:
        for line in inf:
            line = line.strip()
            #skip empty and header line
            if not line or line[0] in header_tokens:
                continue
            else:
                subcl_filt.append(line.split()[0].split(',')[0])
    return subcl_filt

def filt_bgcs_on_subcl_mot(bgcs, filt, exclude_len_1 = True):
    '''Return same bgcs dict only containing filt motifs

    bgcs: {header:{class:'',motif:[info]}}
    filt: list of str, motif names to filter on
    exclude_len_1: bool, indicating of matches of len 1 should be skipped
    '''
    new_bgcs = {}
    cls = 'class'
    for bgc, items in bgcs.items():
        new_items = {}
        for item in items:
            if item == cls:
                new_items[cls] = items[cls]
            else:
                if item in filt:
                    if exclude_len_1:
                        #unpack gene, prob values to list of (gene, prob)
                        g = [(gp.split(':')[0], float(gp.split(':')[1])) for \
                            gp in items[item][-1].split(',')]
                        if not len(g) > 1 and not round(g[0][1]) > 1:
                            continue
                    new_items[item] = items[item]
        if not len(new_items) == 1:
            new_bgcs[bgc] = new_items
    return new_bgcs

def plot_topics_per_bgc(topics_per_bgc, outname):
    '''Make a barplot of the amount of topics per bgc

    topics_per_bgc: dict/counter object, {n:bgcs_with_n_topics}
    outname: str
    '''
    xs = range(max(topics_per_bgc)+1)
    h = [topics_per_bgc[x] if x in topics_per_bgc else 0 for x in xs]
    plt.close()
    plt.bar(xs, h)
    plt.xlabel('Number of annotated topics per BGC')
    plt.ylabel('Occurrence')
    plt.title('Annotated topics per BGC')
    plt.savefig(outname)
    plt.close()


if __name__ == '__main__':
    try:
        bgc_file = argv[1]
        filt_file = argv[2]
    except IndexError:
        print('Incorrect usage. Provide bgc_topics.txt and tsv/csv of '+
            'filtered subcluster motif names, one subcl motif per line.')
        exit()

    bgc_with_subcl = read_txt_to_list(bgc_file)
    print('Amount of intital BGCs:',len(bgc_with_subcl))
    filt_list = read_filter_subcl_mots(filt_file)
    filt_bgcs = filt_bgcs_on_subcl_mot(bgc_with_subcl, filt_list)
    #test:
    #print(filt_bgcs['BGC0001830'])
    #try:
    #    print(filt_bgcs['NC_020410.1.cluster006'])
    #except KeyError:
    #    print('NC_020410.1.cluster006 not found')
    print(len(filt_bgcs) ,'BGCs with one or more sub-cluster motifs from',\
        filt_file)
    num_mibigs = [bgc for bgc in filt_bgcs if bgc.startswith('BGC0')]
    print(len(num_mibigs) ,'MiBIG BGCs with one or more sub-cluster motifs from',\
        filt_file)

    out_plot = bgc_file.strip('.txt') + '-topic_per_bgc_filtered.pdf'
    #-1 for class info
    counts = Counter([len(vals)-1 for vals in filt_bgcs.values()])
    print('\nBarplot at', out_plot)
    plot_topics_per_bgc(counts, out_plot)

    max_count = max(counts.keys())-1
    highest = [bgc for bgc, vals in filt_bgcs.items() if len(vals)-1 >= \
        max_count]
    print('\nBGCs with most annotated topics ({} or more):\n{}'.format(\
        max_count, '\n'.join(highest)))

    other_cls = 'other'
    cls = 'class'
    tot_other = [bgc for bgc, vals in filt_bgcs.items() if vals[cls] ==\
        other_cls]
    other_mibig = [bgc for bgc in tot_other if bgc.startswith('BGC0')]
    print('\n{} BGCs with other class annotation'.format(len(tot_other)))
    print('{} MiBIGs with other class annotation'.format(len(other_mibig)))
    tot_other = [bgc for bgc, vals in filt_bgcs.items() if vals[cls] ==\
        other_cls and bgc[0] in 'NJ']
    #i think this includes all antismashdb BGCs
    print('{} BGCs with other class annotation whose name start with N or J'\
        .format(len(tot_other)-len(other_mibig)))