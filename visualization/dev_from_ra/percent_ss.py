#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate % secondary structure.''')

parser.add_argument('df', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')


def percent_HLS(ss_seq):
	'''Match helix loop and sheet states and return percentages
	'''
	ss = {'G':'H',
          'H':'H',
          'I':'H',
          'E':'S',
          'B':'S',
          'T':'L',
          'S':'L',
          ' ':'L' #COILS
         }

    states = {'H':0, 'L':0, 'S':0}
    counts = Counter(ss_seq)
    for key in ss:
    	states[ss[key]]+=counts[key]

    return states['H']/len(ss_seq), states['L']/len(ss_seq), states['S']/len(ss_seq)


def parse_ss(df, aln_type):
	'''Extract ss and return % H, L, S
	'''
	
	ss_seq1 = [*df['ss1'+aln_type]]
	ss_seq2 = [*df['ss2'+aln_type]]

	#Save helices loops and sheets
	H1 = []
	L1 = []
	S1 = []
	H2 = []
	L2 = []
	S2 = []

	for i in range(len(ss_seq1)):
		#Sequence 1
		ss1 = ss_seq1[i]
		H,L,S = percent_HLS(ss1)
		H1.append(H)
		L1.append(L)
		S1.append(S)
		#Sequence 2
		ss2 = ss_seq2[i]
		H,L,S = percent_HLS(ss2)
		H2.append(H)
		L2.append(L)
		S2.append(S)

	df['H1'+aln_type] = H1
	df['L1'+aln_type] = L1
	df['S1'+aln_type] = S1

	df['H2'+aln_type] = H2
	df['L2'+aln_type] = L2
	df['S2'+aln_type] = S2

	return df


#MAIN
args = parser.parse_args()
df = pd.read_csv(args.df[0])
df = parse_ss(df, '_straln')
