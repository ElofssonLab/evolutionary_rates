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

parser.add_argument('--get_one', nargs=1, type= bool,
default=sys.stdin, help = 'Get one pair from each H-group (1) or all (0).')

def compactness(mldists, scores):
    '''Assess compactness of topology'''

    mldists = mldists/9 #Normalize (lddt is already 0-1)
    X = np.array([mldists,scores])
    Y = pdist(X.T, 'euclidean')
    return Y

def plot_class_distr(df, outdir):
    classes = [1.,2.,3.,4.]
    names = ['Alpha','Beta','AlphaBeta','FewSS']
    for C in classes:
        class_df = df[df['C._x'] == C]
        sns.distplot(class_df['lddt_scores_straln'], label = C)

    plt.legend()
    #plt.show()

    f = open(outdir+'pvals_per_class.txt', 'w')
    f.write('''P-values calculated by two sided t-tests between straln lddt scores wihtin classes using points
            above or equal to 6 in ML AA20 distance, but below or equal to 8.9\n''')

    #Write number of points
    for i in range(0,4):
        Ci =df[df['C._x'] == classes[i]]
        f.write('Number of points for class '+names[i]+': '+str(len(Ci))+'\n')

    #Calculate p-values and write to file
    for i in range(0,3):
        C1 =df[df['C._x'] == classes[i]]
        f.write('\n'+str(classes[i])+'\t'*int(i+2))
        for j in range(i+1,4):
            C2 =df[df['C._x'] == classes[j]]
            statistic, pvalue = stats.ttest_ind([*C1['lddt_scores_straln']], [*C2['lddt_scores_straln']], equal_var = False)
            f.write(str(np.round(pvalue,3))+'\t')
    f.close() #Close file
    return None
#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
top_metrics = pd.read_csv(args.top_metrics[0])
get_one = args.get_one[0]
pdb.set_trace()
cardinality = '_AA20'
#set random seed
np.random.seed(42)
if get_one == True:
    #get one pair per H-group from hgroupdf
    groups = [*Counter(hgroupdf['H_group']).keys()]
    one_pair_df = pd.DataFrame(columns = hgroupdf.columns)
    for g in groups:
        partial_df = hgroupdf[hgroupdf['H_group']==g]
        i = np.random.randint(len(partial_df), size = 1)
        start =  partial_df.index[0]
        selection = partial_df.loc[start+i]
        one_pair_df = one_pair_df.append(selection)
    hgroupdf = one_pair_df

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
partial_df = catdf[catdf['MLAAdist_straln']>=6]
partial_df = partial_df[partial_df['MLAAdist_straln']<=8.9]
#plot_class_distr(partial_df, outdir)


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
    if len(df) <1: #If no points in df
        continue
    mldists = np.asarray(df['MLAAdist_straln'])
    scores = np.asarray(df['lddt_scores_straln'])
    #y = compactness(mldists,scores) #assess spread in topology
    #Y.append(y[0])

    if pval <t:
        if av > 0:
            lab = 'sig+'
        else:
            lab = 'sig-'
    else:
        lab = 'nonsig'


    plt.scatter(mldists, scores, c = colors[lab], s = 3)

plt.legend()
plt.xlabel('ML AA20 distance')
plt.ylabel('lddt_scores_straln')
plt.show()
