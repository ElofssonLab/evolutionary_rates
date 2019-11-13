#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd
import glob
import subprocess
from collections import Counter
import os
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that calculates reduced cardinalities for a df, makes phylip files and runs treepuzzle
                                                to calculate ML AA distances.''')

parser.add_argument('--df', nargs=1, type= str,
                  default=sys.stdin, help = 'path to input dataframe.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'path output directory.')
parser.add_argument('--puzzle', nargs=1, type= str,
                default=sys.stdin, help = '''path to tree-puzzle.''')

def encode2(seq):
    '''Encode aa sequence to H/P
    '''
    AA2 = {
    'L':'H',
    'I':'H',
    'V':'H',
    'F':'H',
    'M':'H',
    'W':'H',
    'C':'H',
    'P':'H',
    'A':'H',
    'K':'P',
    'R':'P',
    'D':'P',
    'E':'P',
    'Y':'P',
    'H':'P',
    'T':'P',
    'S':'P',
    'Q':'P',
    'N':'P',
    'G':'P',
    '-':'-',
    'X':'X'
    }

    aa2seq = ''
    for i in seq:
        aa2seq+=AA2[i]

    return aa2seq

def encode3(seq):
    '''Encode aa sequence to H/P/G
    '''
    AA3 = {
    'L':'H',
    'I':'H',
    'V':'H',
    'F':'H',
    'M':'H',
    'W':'H',
    'C':'H',
    'A':'H',
    'K':'P',
    'R':'P',
    'D':'P',
    'E':'P',
    'Y':'P',
    'H':'P',
    'T':'P',
    'S':'P',
    'Q':'P',
    'N':'P',
    'P':'G',
    'G':'G',
    '-':'-',
    'X':'X'
    }

    aa3seq = ''
    for i in seq:
        aa3seq+=AA3[i]

    return aa3seq

def encode6(seq):
    '''Encode aa sequence to K/D/Y/T/C
    '''
    AA6 = {
    'K':'K',
    'R':'K',
    'D':'D',
    'E':'D',
    'Y':'Y',
    'W':'Y',
    'F':'Y',
    'H':'Y',
    'T':'T',
    'S':'T',
    'Q':'T',
    'N':'T',
    'C':'C',
    'V':'C',
    'M':'C',
    'L':'C',
    'I':'C',
    'A':'C',
    'P':'P',
    'G':'P',
    '-':'-',
    'X':'X'
    }

    aa6seq = ''
    for i in seq:
        aa6seq+=AA6[i]

    return aa6seq

def reduce_cardinality(sequences, suffix, RC_df):

    #Encode sequence
    AA2 = []
    AA3 = []
    AA6 = []
    for s in sequences:
        AA2.append(encode2(s))
        AA3.append(encode3(s))
        AA6.append(encode6(s))

    #New df with AA2,3,6 alphabet
    RC_df['AA2_'+suffix] = AA2
    RC_df['AA3_'+suffix] = AA3
    RC_df['AA6_'+suffix] = AA6

    return RC_df


def make_phylip(uids, query_aln, template_aln, outdir):
        '''Print phylip format for tree-puzzle calculations
        '''
        #Create text in phylip format
        text = (' 4  ' + str(len(query_aln)) + '\n'
                        + uids[0] + '00|' + query_aln + '\n'
                        + 'copy11111' + '|' + query_aln + '\n'
                        + uids[1] + '00|' + template_aln + '\n'
                        + 'copy22222' + '|' + template_aln + '\n')


        #Define file name
        file_name = uids[0] + '_' + uids[1] + '.phy'
        #Open file and write text to it
        with open(outdir+file_name, "w") as file:
                file.write(text)

        return None

def run_puzzle(indir, puzzle):
    '''Run tree-puzzle and retrieve output
    '''

    for name in glob.glob(indir+"*.phy"): #Use all .phy files
        uid_pairs = name.split('/')[-1].split('.')[0].split('_')
        try:
            p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
            p.communicate(b'y\nn\n')[0]
        except:
            raise IOError(name)

    return None

def parse_puzzle(results_dir, RC_df, suffices, cardinalities):
    '''Parse the newly generated tree-puzzle results
    '''
    #base
    base = results_dir+'/reduced_cardinality/'
    uid1 = [*RC_df['uid1']]
    uid2 = [*RC_df['uid2']]
    groups = [*RC_df['group']]

    for cardinality in cardinalities:
        for suffix in suffices:
            evdists = []
            for i in range(len(uid1)):
                key = uid1[i]+'_'+uid2[i]
                indir = base + cardinality + '/' + suffix + '/'+groups[i]+'/'
                #Open distfile - tree-puzzle results
                dist_file = open(indir + key + '.phy.dist')
                for line in dist_file:
                    line = line.rstrip()
                    line = line.split(" ") #split on double space
                    line = list(filter(None, line)) #Filter away empty strings

                    if len(line)>2:
                        seq_dist = line[-1] #Get ML evolutionary distance between sequences
                        break
                dist_file.close()
                evdists.append(seq_dist)
            RC_df['MLAAdist_'+cardinality+'_'+suffix] = evdists
    RC_df.to_csv(results_dir+'RC_df.csv')
    return RC_df

def reduce_and_run(complete_df, results_dir, puzzle):
    #Reduce cardinality, make phylip files and run tree-puzzle
    RC_df = pd.DataFrame()
    #Encode sequences with reduced cardinality
    suffices = ['1_seqaln', '2_seqaln', '1_straln', '2_straln']
    for suffix in suffices:
        sequences = complete_df['seq'+suffix]
        RC_df = reduce_cardinality(sequences, suffix, RC_df)


    RC_df['uid1'] = complete_df['uid1']
    RC_df['uid2'] = complete_df['uid2']
    RC_df['group'] = complete_df['group']

    RC_df.to_csv(results_dir+'RC_df.csv') #Save RC df
    #Write .phy files for reduced cardinality representations of alignments
    suffices = ['seqaln', 'straln']
    cardinalities = ['AA2', 'AA3', 'AA6']
    u_groups = [*Counter(RC_df['group']).keys()]
    os.mkdir(results_dir+'/reduced_cardinality/')
    for cardinality in cardinalities:
        os.mkdir(results_dir+'/reduced_cardinality/'+cardinality)
        for suffix in suffices:
            os.mkdir(results_dir+'/reduced_cardinality/'+cardinality+'/'+suffix)
            for group in u_groups:
                outdir = results_dir+'/reduced_cardinality/'+cardinality+'/'+suffix+'/'+group+'/'
                if not os.path.isdir(outdir): #check if dir exists, otherwise make it
                    os.mkdir(outdir)

                df = RC_df[RC_df['group'] == group]

                seq1 = [*df[cardinality+'_1_'+suffix]]
                seq2 = [*df[cardinality+'_2_'+suffix]]
                uid1 = [*df['uid1']]
                uid2 = [*df['uid2']]
                for i in range(len(seq1)):
                    make_phylip([uid1[i], uid2[i]], seq1[i], seq2[i], outdir)


    #Run tree-puzzle on reduced cardinality phylip files
    for cardinality in cardinalities:
        for suffix in suffices:
             print(cardinality,suffix)
             dirs = glob.glob(results_dir+'/reduced_cardinality/'+cardinality+'/'+suffix+'/*')
             for d in dirs:
                 run_puzzle(d+'/', puzzle)

    #Parse tree-puzzle results
    RC_df = parse_puzzle(results_dir, RC_df, suffices, cardinalities)
    #Add columns to complete_df
    columns = ['AA2_1_seqaln', 'AA3_1_seqaln', 'AA6_1_seqaln', 'AA2_2_seqaln', 'AA3_2_seqaln',
           'AA6_2_seqaln', 'AA2_1_straln', 'AA3_1_straln', 'AA6_1_straln', 'AA2_2_straln',
            'AA3_2_straln', 'AA6_2_straln', 'MLAAdist_AA2_straln', 'MLAAdist_AA3_straln',
           'MLAAdist_AA6_straln', 'MLAAdist_AA2_seqaln', 'MLAAdist_AA3_seqaln',
           'MLAAdist_AA6_seqaln']

    for column in columns:
        complete_df[column] = RC_df[column]

    complete_df.to_csv(results_dir+'complete_df_allcards.csv')

#####MAIN#####
args = parser.parse_args()

df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
puzzle = args.puzzle[0]

#Reduce cardinalities
reduce_and_run(df, outdir, puzzle)
