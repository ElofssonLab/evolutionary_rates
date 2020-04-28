#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
from collections import Counter
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Read in ligands for all pdb entries and format matching/non-matching sets based on the ligands of each pair.''')

parser.add_argument('--ligands_per_pdb', nargs=1, type= str, default=sys.stdin, help = 'Path to ligands per pdb id.')

parser.add_argument('--topdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--av_df', nargs=1, type= str, default=sys.stdin, help = 'Path to df with running averages.')

parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

###FUNCTIONS###
def parse_ligands(ligand_file):
    '''Get all ligands for each pdb id into a dict
    '''
    ligands_per_pdb = {} #Connection btw pdb id and its ligands
    with open(ligand_file, 'r') as file:
        for line in file:
            line = "".join(line.split(" ")) #Remove space
            line = line.rstrip() #lagging newlines
            line = line.split(':')
            pdb_id = line[0]
            ligands = line[1].split(';')[:-1]
            ligands_per_pdb[pdb_id]=ligands

    return ligands_per_pdb

def match_ligands(catdf, ligands_per_pdb):
    '''Match ligands for each uid in all pairs in df
    '''
    same_ligand = [] #Indicates if the pdb ids in each pairs have the same ligands(1) or not (0)
    max_shared_ligands_in_pair = []
    for i in range(len(catdf)):
        row = catdf.iloc[i]
        try:
            lig1 = ligands_per_pdb[row['uid1'][:4]]
            lig2 = ligands_per_pdb[row['uid2'][:4]]
        except:
            same_ligand.append('NA')
            max_shared_ligands_in_pair.append('NA')
            continue
        #Check if they share
        num_shared = 0
        for lig in lig1:
            if lig in lig2:
                num_shared +=1
        same_ligand.append(num_shared)
        #Append min num (max possible shared)
        max_shared_ligands_in_pair.append(min([len(lig1),len(lig2)]))

    #Add to catdf
    catdf['shared_ligands'] = same_ligand
    catdf['max_shared_ligands_in_pair']=max_shared_ligands_in_pair
    return catdf


#####MAIN#####
args = parser.parse_args()
ligands_per_pdb = parse_ligands(args.ligands_per_pdb[0])
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
av_df = pd.read_csv(args.av_df[0])
outdir = args.outdir[0]

cardinality = '_AA20'

#Get topology from hgroupdf
tops = []
hgroups = [*hgroupdf['group']]
for hg in hgroups:
    hg = hg.split('.')
    tops.append(hg[0]+'.'+hg[1]+'.'+hg[2])

hgroupdf['C.A.T.'] = tops
#rename col
hgroupdf = hgroupdf.rename(columns={'group':'H_group'})
hgroupdf = hgroupdf.rename(columns={'C.A.T.':'group'})

#rename TMscore cols
hgroupdf = hgroupdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})
topdf = topdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})

catdf = pd.concat([topdf, hgroupdf])
#Rename class column
catdf = catdf.rename(columns={'C._x':'Class'})

#The ones should actually be zeros
catdf['RCO1']=catdf['RCO1'].replace([1], 0)
catdf['RCO2']=catdf['RCO2'].replace([1], 0)

#Match ligands in each pair
catdf = match_ligands(catdf, ligands_per_pdb)
