#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import subprocess
import numpy as np
import sys
import os
import argparse
import glob
import pdb

#custom imports
sys.path.insert(0, '/home/pbryant/evolutionary_rates/')
print(sys.path)
from conversions import make_phylip


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle on domains from SCOP.''')

parser.add_argument('--df', nargs=1, type= str, default=sys.stdin, help = 'path to df.')
parser.add_argument('--indir', nargs=1, type= str, default=sys.stdin, help = 'Path to input directory.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('--puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')
parser.add_argument('--TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')


#####FUNCTIONS#####
def align(df, indir, outdir, TMalign):
        '''Run TMalign on .pdb files
        '''
        measures = {} #Save RMSD to add with MLAA distance from tree-puzzle

        for i in range(len(df)):
            row = df.iloc[i]

            #Run TMalign and extract scores
            str1 = indir+row['uid1']+'.ent'
            str2 = indir+row['uid2']+'.ent'
            try:
                tmalign_out = subprocess.check_output([TMalign, str1 , str2]) #Performs optimal structural alignment
            except:
                continue
            (tm_aligned_len, rmsd, tmscores, tm_identity, chain_lens, tm_sequences)= parse_tm(tmalign_out)
            measures[row['uid1']+'_'+row['uid2']] = [rmsd, tmscores[0], tmscores[1]]

            #Write .phy file of alignment
            make_phylip([row['uid1'], row['uid2']], tm_sequences[0], tm_sequences[1], outdir)

        return measures

def parse_tm(tmalign_out):
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

def run_puzzle(outdir, puzzle):
        '''Run tree-puzzle and retrieve output
        '''
        for name in glob.glob(outdir+"*.phy"): #Use all .phy files
                try:
                        p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
                        p.communicate(b'y\nn\n')[0]
                except:
                        raise IOError(name)


        return None

def parse_puzzle(measures, indir):
        '''Parse output from tree-puzzle and write to dict
        '''
        keys = [*measures] #Make list of keys in dict
        for key in keys:
                rmsd, tmscore1, tmscore2 = measures[key] #Get rmsd
                try:
                        dist_file = open(indir + key + '.phy.dist')
                except:
                        uids = [key[0:7],key[8:]]
                        dist_file = open(indir + uids[1] + '_' + uids[0] + '.phy.dist')
                        measures.pop(key)
                        #change key to match other file names
                        key = uids[1] + '_' + uids[0]
                for line in dist_file:
                        line = line.rstrip()
                        line = line.split(" ") #split on double space
                        line = list(filter(None, line)) #Filter away empty strings

                        if len(line)>2:
                                seq_dist = line[-1] #Get ML evolutionary distance between sequences
                                measures[key] = [rmsd, tmscore1, tmscore2, seq_dist]
                                break
                dist_file.close()

        return measures

def print_tsv(measures, hgroup, outdir):
        '''Print measures in tsv to file
        '''
        with open(outdir+hgroup+'_str.csv', 'w') as file:
                file.write('uid1,uid2,MLAAdist,RMSD,TMscore_high,TMscore_low\n')
                for key in measures:
                    uid1 = key[0:7]
                    uid2 = key[8:]
                    rmsd, tmscore1, tmscore2, seq_dist = measures[key]
                    high_score = max(float(tmscore1), float(tmscore2))
                    low_score = min(float(tmscore1), float(tmscore2))
                    file.write(uid1+','+uid2+','+seq_dist+','+rmsd+','+str(high_score)+','+str(low_score)+'\n')

        return None

#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
indir = args.indir[0]
outdir = args.outdir[0]
puzzle = args.puzzle[0]
TMalign = args.TMalign[0]

#Run TMalign
measures = align(df, indir, outdir, TMalign)
print('Made alignments with TMalign')
#Run puzzle
run_puzzle(outdir, puzzle)
print('Run tree-puzzle')

measures = parse_puzzle(measures, outdir)
print_tsv(measures, '2009_SCOP', outdir)
