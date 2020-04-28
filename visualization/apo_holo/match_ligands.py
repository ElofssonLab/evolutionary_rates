#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np

import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Read in ligands for all pdb entries and format matching/non-matching sets based on the ligands of each pair.''')

parser.add_argument('--ligands_per_pdb', nargs=1, type= str, default=sys.stdin, help = 'Path to ligands per pdb id.')

parser.add_argument('--topdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')



#####MAIN#####
args = parser.parse_args()
ligands_per_pdb = parse_ligands(args.ligands_per_pdb[0])
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
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

