#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that collects results and merges them
                                                into a unified datframe. A selection on how much of
                                                the shortest sequence in each pair is made according
                                                to a threshold.''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to input directory.')
parser.add_argument('threshold', nargs=1, type= float,
                  default=sys.stdin, help = 'Threshold for percent aligned of the shortest sequence in each pair.')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'path output directory.')


def tsv_to_df(folder_path):
    '''Get results from tsv files into 
    '''
    #Get results from sequence alignments
    all_files = glob.glob(results_path+'*/TMscore/*seq.tsv')     # advisable to use os.path.join as this makes concatenation OS independent
    df_from_each_file = [pd.read_csv(f, sep='\t') for f in all_files]
    for dataframe, filename in zip(df_from_each_file, all_files):
    hgroup = filename.split('/')[-3]
    dataframe['H_group'] = hgroup
    dataframe['C.'] = hgroup.split('.')[0]+'.'
    dataframe['C.A.'] = hgroup.split('.')[0]+'.'+hgroup.split('.')[1]
    from_seq_df = pd.concat(df_from_each_file, ignore_index=True)
    from_seq_df.to_csv(results_path+'/from_seq_df.csv')

def create_df(results_path):
    '''Gets results and parses them into a unified dataframe
    '''

    #Get results from sequence alignments
    all_files = glob.glob(results_path+'*/TMscore/*seq.tsv')     # advisable to use os.path.join as this makes concatenation OS independent
    df_from_each_file = [pd.read_csv(f, sep='\t') for f in all_files]
    for dataframe, filename in zip(df_from_each_file, all_files):
    hgroup = filename.split('/')[-3]
    dataframe['H_group'] = hgroup
    dataframe['C.'] = hgroup.split('.')[0]+'.'
    dataframe['C.A.'] = hgroup.split('.')[0]+'.'+hgroup.split('.')[1]
    from_seq_df = pd.concat(df_from_each_file, ignore_index=True)
    from_seq_df.to_csv(results_path+'/from_seq_df.csv')




#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
t = args.threshold[0]
outdir = args.outdir[0]
