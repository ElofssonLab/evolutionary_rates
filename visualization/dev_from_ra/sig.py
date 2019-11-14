#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
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

parser.add_argument('--get_one', nargs=1, type= int,
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
            f.write(str(pvalue)+'\t')
    f.close() #Close file
    return None


def plot_rco(catdf, aln_type, cmap_name):
    #Plot RCO
    matplotlib.rcParams.update({'font.size': 22})

    fig = plt.figure(figsize=(12,10)) #set figsize
    cmap = cm.get_cmap(cmap_name, 10) #twilight_shifted was wuite nice - donät see diff btw low and high though
    step = 0.1
    for t in np.arange(step, 1+step, step):
        partial_rco1 = catdf[catdf['RCO1']<=t]
        partial_rco1 = partial_rco1[partial_rco1['RCO1']>t-step] #RCO1 in interval

        partial_rco2 = catdf[catdf['RCO2']<t]
        partial_rco2 = partial_rco2[partial_rco2['RCO2']>=t-step] #RCO2 in interval
        #concat
        cat_rco = pd.concat([partial_rco1, partial_rco2])
        cat_rco = cat_rco.drop_duplicates() #Drop duplicates
        plt.scatter(cat_rco['MLAAdist'+aln_type], cat_rco['lddt_scores'+aln_type], color = cmap(t), label = str(t), s=2, alpha = 1)

    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    #cbar.set_label('Relative contact order', rotation=270)
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score')
    plt.show()
    fig.savefig(outdir+aln_type+'_RCO.svg', format = 'svg')
    return None

def outliers(catdf):
    catdf = catdf[catdf['lddt_scores_straln']>0.9]
    catdf = catdf[catdf['MLAAdist_straln']>2]
    uids = [*catdf['uid1']]+[*catdf['uid2']]
    uids = Counter(uids)


    pdb.set_trace()

def plot_sig(catdf, top_metrics):
    '''Plot pairs py siginificance
    '''
    Y = []
    t = 0.05/len(top_metrics)
    sig_rco = 0 #Save assessment if sig and rco >0.9 or <0.1
    colors = {'negative':'g', 'non-significant':'k', 'positive':'b'}
    matplotlib.rcParams.update({'font.size': 22})
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
            #Assess RCO
            rco1 = np.array(df['RCO1'])
            rco2 = np.array(df['RCO2'])
            rco = np.concatenate([rco1, rco2])

            if len(np.where(rco>0.9)[0]) > 0 or len(np.where(rco<0.1)[0])>0:
                sig_rco = 1
            else:
                sig_rco = 0

            if sig_rco == 1:
                lab = 'positive'
                alpha = 0.8
            else:
                lab = 'negative'
                alpha = 0.8

            # if av > 0:
            #     lab = 'positive'
            #     alpha = 1
            # else:
            #     lab = 'negative'
            #     alpha = 1
        else:
            lab = 'non-significant'
            alpha = 0.8


        plt.scatter(mldists, scores, c = colors[lab], s = 3, alpha = alpha)

    #plt.legend(('negative', 'non-significant', 'positive'),
               #shadow=False, loc=(0.48, 0.75), fontsize=22, markerscale=5., scatterpoints=1)
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score')
    plt.show()
    fig.savefig(outdir+'lddt_straln_rco_and_sig.svg', format = 'svg')

    return None
#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
top_metrics = pd.read_csv(args.top_metrics[0])
get_one = bool(args.get_one[0])

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
hgroups = [*hgroupdf['group']]
for hg in hgroups:
    hg = hg.split('.')
    tops.append(hg[0]+'.'+hg[1]+'.'+hg[2])

hgroupdf['C.A.T.'] = tops

#rename col
hgroupdf = hgroupdf.rename(columns={'group':'H_group'})
hgroupdf = hgroupdf.rename(columns={'C.A.T.':'group'})
catdf = pd.concat([topdf, hgroupdf])

#The ones should actually be zeros
catdf['RCO1']=catdf['RCO1'].replace([1], 0)
catdf['RCO2']=catdf['RCO2'].replace([1], 0)


#Plot by RCO
plot_rco(catdf, '_straln', 'viridis')
pdb.set_trace()
#sig_and_rco(catdf)
partial_df = catdf[catdf['MLAAdist_straln']>=6]
partial_df = partial_df[partial_df['MLAAdist_straln']<=8.9]
#plot_class_distr(partial_df, outdir)

plot_sig(catdf, top_metrics)
