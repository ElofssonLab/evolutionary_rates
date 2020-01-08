#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import subprocess
import numpy as np
import sys
import os
import argparse
import pdb

#custom imports
sys.path.insert(1, '/home/pbryant/evolutionary_rates')
from conversions import make_phylip
from run_tmalign_treepuzzle_ind import run_puzzle

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb files for the old alignments
                                                for running the new flow on them.''')

parser.add_argument('--df', nargs=1, type= str, default=sys.stdin, help = 'path to df.')
parser.add_argument('--indir', nargs=1, type= str, default=sys.stdin, help = 'Path to input directory.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('--TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')


#####FUNCTIONS#####
def run_TMalign(df, indir, outdir, TMalign):
	'''Run TMalign on .pdb files
	'''

    for i in range(len(df)):
        row = df.iloc[i]

        #Run TMalign and extract scores
        str1 = indir+row['uid1']+'.ent'
        str2 = indir+row['uid2']+'.ent'

        tmalign_out = subprocess.check_output([TMalign, str1 , str2]) #Performs optimal structural alignment
        (tm_aligned_len, rmsd, tmscores, tm_identity, chain_lens, tm_sequences)= parse_tm(tmalign_out)
        measures[uid1+'_'+uid2] = [rmsd, tmscores[0], tmscores[1]]

        #Write .phy file of alignment
        make_phylip(uids, tm_sequences[0], tm_sequences[1], outdir)

ef parse_tm(tmalign_out):
	'''A function that parses TMalign output.
	'''

	tmalign_out = tmalign_out.decode("utf-8")
	tmalign_out = tmalign_out.split('\n')
	tmscores = [] #Save TMscores
	for i in range(0, len(tmalign_out)): #Step through all items in list

		if 'Aligned length' and 'RMSD' and 'Seq_ID' in tmalign_out[i]:
			row = tmalign_out[i].split(',')
			aligned_len = row[0].split('=')[1].lstrip()
			rmsd = row[1].split('=')[1].lstrip()
			identity = row[2].split('=')[2].lstrip()

		if 'Length of Chain_1:' in tmalign_out[i]:
			len_1 = tmalign_out[i].split(':')[1].split()[0]

		if 'Length of Chain_2:' in tmalign_out[i]:
                        len_2 = tmalign_out[i].split(':')[1].split()[0]
		if 'TM-score=' in tmalign_out[i]:
			tmscores.append(tmalign_out[i].split('(')[0].split('=')[1].strip())

	#Get per residue sequence alignments from structural alignment
	sequences = [tmalign_out[-5], tmalign_out[-3]]

	chain_lens = [int(len_1), int(len_2)]


	return(aligned_len, rmsd, tmscores, identity, chain_lens, sequences)
    
#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
