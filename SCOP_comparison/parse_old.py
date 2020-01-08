#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
from collections import Counter
import subprocess
import numpy as np
import sys
import os
import argparse
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses the 2009 alignments and returns
                                                the SCOP identifiers as well as alignments for tree-puzzle.''')

parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'path to infile.')

#####FUNCTIONS#####
def parse(infile):
    '''A function that parses the alignments from the 2009 dataset
    Illergård, K., Ardell, D. H. & Elofsson, A. Structure is three to ten times more conserved than sequence-
    A study of structural response in protein cores.
    Proteins: Structure, Function, and Bioinformatics vol. 77 499–508 (2009).
    '''

    #Save information about ids and alignments
    uid1 = []
    uid2 = []
    alignments = [[], []] #refers to alignments of all AA20 (all standard amino acids)

    with open(infile, 'r') as f:
        get_next = False #Keep track of alignments fetched
        for line in f:
            line = line.rstrip()
            if '= vs =' in line:
                line = line.split('=')
                uid1.append(line[13].strip()) #some contain '_' in end, this should apparently be there
                uid2.append(line[23].strip())
            if '>Ali1_AA' == line:
                get_next = True
                index = 0
                continue

            if '>Ali2_AA' == line:
                get_next = True
                index = 1
                continue

            if get_next == True:
                alignments[index].append(line)
                get_next = False


    df = pd.DataFrame()
    df['uid1'] = uid1
    df['uid2'] = uid2
    df['aln1'] = alignments[0]
    df['aln2'] = alignments[1]
    df.to_csv('alignment_df.csv')

    uids = uid1+uid2
    return df, uids

def get_pdbs(uids):
    '''A script for getting pdb files
    '''
    f = open('failed_ids.txt', 'a+')
    base='https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-1.75/'
    for uid in uids:
        if not os.path.isfile(uid+'.ent'):
            try:
                outp = subprocess.check_output('wget '+base+uid[2:4]+'/'+uid+'.ent', shell=True)
            except:
                f.write(uid+'\n')

#####MAIN#####
args = parser.parse_args()
infile = args.infile[0]

df, uids = parse(infile)
#Get unique uids
uids = [*Counter(uids).keys()]
with open('uids.txt', 'w') as f:
    for id in uids:
        f.write(id+'\n')
print('Parsed')

get_pdbs(uids)
