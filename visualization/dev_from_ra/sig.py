#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy import stats
from scipy.spatial.distance import pdist

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that .''')

parser.add_argument('--topdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')

parser.add_argument('--top_metrics', nargs=1, type= str,
default=sys.stdin, help = 'Dataframe with topology grouping metrics.')

def compactness(mldists, scores):
    '''Assess compactness of topology'''

    mldists = mldists/9 #Normalize (lddt is already 0-1)
    X = np.array([mldists,scores])
    Y = pdist(X.T, 'euclidean')
    return Y

#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
top_metrics = pd.read_csv(args.top_metrics[0])

cardinality = '_AA20'

#Get topology from hgroupdf
tops = []
hgroups = [*hgroupdf['H_group']]
for hg in hgroups:
    hg = hg.split('.')
    tops.append(hg[0]+'.'+hg[1]+'.'+hg[2])

hgroupdf['C.A.T.'] = tops
#rename col
hgroupdf = hgroupdf.rename(columns={'C.A.T.':'group'})
catdf = pd.concat([topdf, hgroupdf])
Y = []
t = 0.05/len(top_metrics)
colors = {'sig+':'r', 'nonsig':'k', 'sig-':'b'}
fig = plt.figure(figsize=(10,10)) #set figsize
for i in range(len(top_metrics)):
    row = top_metrics.iloc[i]
    top = row['Topology']
    pval = row['lddt_scores_straln_pval']
    av = row['lddt_scores_straln_av_dev']
    df = catdf[catdf['group']==top]
    mldists = np.asarray(df['MLAAdist_straln'])
    scores = np.asarray(df['lddt_scores_straln'])
    y = compactness(mldists,scores)
    Y.append(y[0])

    if pval <t:
        if av > 0:
            lab = 'sig+'
        else:
            lab = 'sig-'
    else:
        lab = 'nonsig'


    plt.scatter(mldists[0:10], scores[0:10], c = colors[lab], s = 3)

plt.legend()
plt.xlabel('ML AA20 distance')
plt.ylabel('lddt_scores_straln')
plt.show()
