#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to run LDA multiple times to find good parameters for applying LDA
to identify modules.
"""

import subprocess

if __name__ == '__main__':
    command = 'python3 lda.py -i /mnt/scratch/louwe015/output_antismashdb_crusemann_mibig_biosynt_subset_19-6/antismashdb_crusemann_mibig_filtered_clusterfile.csv -m /mnt/scratch/louwe015/output_antismashdb_crusemann_mibig_biosynt_subset_19-6/antismashdb_crusemann_mibig_filtered_modules.txt -o /mnt/scratch/louwe015/lda_antismashdb_crusemann_mibig_biosynt_subset_{}_{}27-6_test -C {} -c 10 -t {} {}-f 0.95 -n 15 --classes /mnt/scratch/louwe015/2018-07-13_16-23-26_hybrids_global_crusemann_bgcs_globalmode_mix_mibig_tsvs/big-scape_classes.txt --known_subclusters /mnt/scratch/louwe015/subclusters_subclusterblast_domains_synt_subset.txt --min_genes 2 -I 10'
    topic_range= [350,500,600,750,1000]
    for i in topic_range:
        #without -a
        command_without_a = command.format(i,'','6000',i,'')
        print(command_without_a)
        try:
            subprocess.check_call(command_without_a, shell=True)
        except subprocess.CalledProcessError:
            print(command_without_a)
            print('all empty topics?')
        
        command_with_a = command.format(i,'x10_','60000',i,'-a 10 ')
        print(command_with_a)
        try:
            subprocess.check_call(command_with_a, shell=True)
        except subprocess.CalledProcessError:
            print(command_with_a)
            print('all empty topics?')
