#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to run LDA multiple times to find good parameters for applying LDA
to identify modules.
"""

import subprocess

if __name__ == '__main__':
    command = 'python3 ~/thesis/scripts-thesis/presto_top/presto_top.py -i ipresto_output_build_models_all_sponge_gbks_biosynt_15-1/all_sponge_gbks_clusterfile_filtered_clusterfile.csv -o ipresto_output_build_models_all_sponge_gbks_biosynt_15-1/ipresto_out_100t_1000chnk -t 100 -C 1000 -I 2000 -c 5 --classes ipresto_output_all_sponge_gbks_10-1/all_classes.txt --known_subclusters /mnt/scratch/louwe015/subcluster_data/subclusterblast_data/subclusters_subclusterblast_domains_synt_subset.txt | tee log_ipresto_output_build_models_all_sponge_gbks_biosynt_15-1_ipresto_out_100t_1000chnk.txt'
    topic_range= [50,75,100,150,200]
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
