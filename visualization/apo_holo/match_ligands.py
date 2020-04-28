#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Read in ligands for all pdb entries and format matching/non-matching sets based on the ligands of each pair.''')

parser.add_argument('--ligands_per_pdb', nargs=1, type= str, default=sys.stdin, help = 'Path to ligands per pdb id.')

parser.add_argument('--topdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str, default=sys.stdin, help = 'path to df.')

parser.add_argument('--avdf', nargs=1, type= str, default=sys.stdin, help = 'Path to df with running averages.')

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
            same_ligand.append(0)#One of the structures does not have any ligands
            max_shared_ligands_in_pair.append(0)
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

def compare_ligand_effect(catdf, avdf):
    '''Compare the effects of having shared bound ligands or not
    '''

    share_ligands = catdf[catdf['shared_ligands']>=1]
    diff_ligands = catdf[catdf['shared_ligands']==0]    

    #Calculate RA for the pairs with shared and different ligands
    start = 0
    end = 6
    step = 0.1
    cardinality = ''
    score = 'lddt_scores'
    aln_type = '_straln'
    #Save averages
    share_av = []
    diff_av = []
    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_share= share_ligands[share_ligands['MLAAdist'+cardinality+aln_type]<=j]
            below_diff= diff_ligands[diff_ligands['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_share= share_ligands[share_ligands['MLAAdist'+cardinality+aln_type]<j]
            below_diff= diff_ligands[diff_ligands['MLAAdist'+cardinality+aln_type]<j]

        below_share = below_share[below_share['MLAAdist'+cardinality+aln_type]>=j-step]
        below_diff = below_diff[below_diff['MLAAdist'+cardinality+aln_type]>=j-step]

        cut_share = np.asarray(below_share[score+aln_type])
        cut_diff = np.asarray(below_diff[score+aln_type])
        #Get averages in interval
        share_av.append(np.average(cut_share))
        diff_av.append(np.average(cut_diff))
   
    #plot
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54))
    x = avdf['ML  distance'].loc[0:59]
    #Plot total av
    ax.plot(x, avdf[score+aln_type].loc[0:59], color = 'darkblue', linewidth = 1, label = 'All ligands')
    #Plot share av
    ax.plot(x, share_av, color = 'g', linewidth = 1, label = 'Shared ligands')
    sns.kdeplot(share_ligands['MLAAdist'+aln_type], share_ligands[score+aln_type], shade = False, cmap = 'Greens')
    #Plot diff av
    sns.kdeplot(diff_ligands['MLAAdist'+aln_type], diff_ligands[score+aln_type], shade = False, cmap = 'Blues')
    ax.plot(x, diff_av, color = 'cornflowerblue', linewidth = 1, label = 'Different ligands')
    #Format plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('AA20 ED')
    ax.set_ylabel('lDDT score')
    ax.set_xlim([0,6.1])
    ax.set_xticks([0,1,2,3,4,5,6])
    ax.set_ylim([0.2,1])
    fig.legend()
    fig.tight_layout()
    fig.savefig(outdir+'apo_holo.png', format = 'png', dpi = 300)
    plt.close()



    pdb.set_trace()
#####MAIN#####
args = parser.parse_args()
ligands_per_pdb = parse_ligands(args.ligands_per_pdb[0])
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
avdf = pd.read_csv(args.avdf[0])
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
#Compare apo and holo forms
compare_ligand_effect(catdf, avdf)
