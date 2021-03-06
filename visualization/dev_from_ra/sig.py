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

parser.add_argument('--avdf', nargs=1, type= str,
default=sys.stdin, help = 'Dataframe with running averages for different scores for topdf+hgroupdf.')

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

    f = open(outdir+'pvals_per_class_1_per_group.txt', 'w')
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


def plot_rco(catdf, aln_type, cmap_name, type, js, avs):
    #Plot RCO
    matplotlib.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(12,10)) #set figsize
    cmap = cm.get_cmap(cmap_name, 10) #twilight_shifted was wuite nice - donät see diff btw low and high though
    step = 0.1
    for t in np.arange(step, 0.7+step, step):
        partial_rco = catdf[catdf['RCO'+type]<=t]
        partial_rco = partial_rco[partial_rco['RCO'+type]>t-step] #RCO1 in interval
        plt.scatter(partial_rco['MLAAdist'+aln_type],partial_rco['lddt_scores'+aln_type], color = cmap(t), s=2, alpha = 0.8)

    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=np.arange(0,0.8,0.1), label = 'RCO'+type)


    #Plot RA as well
    plt.plot(js, avs, linewidth = 2, c = 'b', label = 'Running average')
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score')
    plt.legend()
    plt.show()
    fig.savefig(outdir+'RCO'+type+'_'+cmap_name+'_lddt_scores'+aln_type+'_RCO.png', format = 'png')
    return None

def RCO_vs_deviation(catdf, avdf, outdir, type):
    '''Plot RCO against the deviation for each point
    '''

    matplotlib.rcParams.update({'font.size': 22})
    #Plot RCO 1 and 2 distibutions
    if type == '1':
        fig = plt.figure(figsize=(10,10)) #set figsize
        sns.distplot(catdf['RCO1'], label = 'RCO1')
        sns.distplot(catdf['RCO2'], label = 'RCO2')
        plt.xlabel('RCO')
        plt.ylabel('Density')
        plt.legend()
        plt.annotate('Pearson R: '+str(np.round(stats.pearsonr(catdf['RCO1'],catdf['RCO2'])[0],2)), xy=(0.4,4))
        plt.annotate('Average RCO: '+str(np.round(np.average(np.array(catdf['RCO1']+catdf['RCO2'])/2),2)), xy=(0.4,3))
        fig.savefig(outdir+'RCOdistributions.svg', format = 'svg')

    plt.close()
    fig = plt.figure(figsize=(10,10)) #set figsize
    #Select RCOs
    deviations = {'low':[], 'high':[], 'middle':[]}
    RCOs = {'low':[], 'high':[], 'middle':[]}
    labels = {'low':'RCO'+type+'<0.1', 'high':'RCO'+type+'>0.5', 'middle':'RCO'+type+' 0.1-0.5'}
    step = 0.1
    for j in np.arange(step,6+step,step):
        tav = avdf[avdf['ML  distance']==np.round(j-0.05,2)]['lddt_scores_straln'].values[0] #total running average
        below_df = catdf[catdf['MLAAdist_straln']<j]
        below_df = below_df[below_df['MLAAdist_straln']>=j-step]
        rcolow = below_df[below_df['RCO'+type]<0.1]
        rcohigh = below_df[below_df['RCO'+type]>0.5]
        rcomiddle = below_df[below_df['RCO'+type]>=0.1]
        rcomiddle = rcomiddle[rcomiddle['RCO'+type]<=0.5]
        RCOs['low'].extend(rcolow['RCO'+type])
        RCOs['high'].extend(rcohigh['RCO'+type])
        RCOs['middle'].extend(rcomiddle['RCO'+type])
        deviations['low'].extend(rcolow['lddt_scores_straln']-tav)
        deviations['high'].extend(rcohigh['lddt_scores_straln']-tav)
        deviations['middle'].extend(rcomiddle['lddt_scores_straln']-tav)

    for key in deviations:
        sns.kdeplot(RCOs[key], deviations[key], shade=True, shade_lowest = False, label = labels[key])
    plt.legend()
    plt.ylabel('Deviation from total running average')
    plt.xlabel('RCO'+type)
    plt.yticks(np.arange(-0.4,0.5,0.1))
    #plt.xticks(np.arange(0,0.7,0.1))
    fig.savefig(outdir+'RCO'+type+'_and_dev_from_ra.svg', format = 'svg')

    #Plot distributions of deviations
    fig = plt.figure(figsize=(10,10)) #set figsize
    sns.distplot(deviations['low'], label = 'RCO'+type+' < 0.1', norm_hist = True)
    sns.distplot(deviations['high'], label = 'RCO'+type+' > 0.5', norm_hist = True)
    sns.distplot(deviations['middle'], label = 'RCO'+type+' 0.1-0.5', norm_hist = True)
    plt.legend()
    plt.xlabel('Deviation from total running average')
    plt.ylabel('Density')
    fig.savefig(outdir+'deviation_distribution_by_RCO'+type+'.svg', format = 'svg')

    f = open(outdir+'pvals_comparing_dev_by_rco'+type+'.txt', 'w')
    f.write('P-values calculated by two sided t-tests between straln lddt scores with rco'+type+' values\n')
    f.write('<0.1, 0.1-0.5 and >0.5 in range 0.6 ML AA20 distance\n')

    keys = ['low', 'high', 'middle']
    f.write('Number of pairs:\n')
    for key in keys:
        f.write(key+': '+str(len(deviations[key]))+'\t')

    for i in range(3):
        key1 = keys[i]
        f.write('\n'+key1+'\t'*int(i+2))
        for j in range(i+1,3):
            key2= keys[j]
            statistic, pvalue = stats.ttest_ind(deviations[key1], deviations[key2], equal_var = False)
            f.write(str(pvalue)+'\t')
    f.close()
    return None

def outliers(catdf, type, outdir):
    if type == 'pos':
        catdf = catdf[catdf['lddt_scores_straln']>0.9]
        catdf = catdf[catdf['MLAAdist_straln']>2]
        uids = [*catdf['uid1']]+[*catdf['uid2']]
        uids = Counter(uids)
    if type == 'neg':
        catdf = catdf[catdf['lddt_scores_straln']<0.45]
        catdf = catdf[catdf['MLAAdist_straln']<2]
        uids = [*catdf['uid1']]+[*catdf['uid2']]
        uids = Counter(uids)

    #Plot RCOs
    matplotlib.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(10,10)) #set figsize
    sns.distplot(catdf['RCO1'], label = 'RCO1', bins = 15)
    sns.distplot(catdf['RCO2'], label = 'RCO2', bins = 15)
    plt.axvline(x=0.29, linestyle='--', linewidth = 2, label = 'Average RCO')
    plt.legend()
    plt.xlabel('RCO')
    plt.ylabel('Density')
    plt.show()

    fig.savefig(outdir+'RCO_distributions_for_'+type+'_outliers.svg', format = 'svg')

def plot_sig(catdf, top_metrics):
    '''Plot pairs by siginificance
    '''
    t = 0.05/len(top_metrics) #Number of groups tested together
    sig_rco = 0 #Save assessment if sig and rco >0.9 or <0.1
    colors = {'negative':'b', 'non-significant':'k', 'positive':'r'}
    alpha = {'negative':1, 'non-significant':0.5, 'positive':1}
    mldists = {'negative':[], 'non-significant':[], 'positive':[]}
    scores = {'negative':[], 'non-significant':[], 'positive':[]}
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

        #y = compactness(mldists,scores) #assess spread in topology
        if pval <t:
            if av > 0:
                lab = 'positive'
            else:
                lab = 'negative'
        else:
            lab = 'non-significant'

        mldists[lab].extend([*df['MLAAdist_straln']])
        scores[lab].extend([*df['lddt_scores_straln']])
    for lab in colors:
        plt.scatter(mldists[lab], scores[lab], c = colors[lab], s = 4, alpha = 0.3, label = lab)

    plt.legend(markerscale=10)
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score')
    plt.show()
    fig.savefig(outdir+'one_pair_lddt_straln_sig.svg', format = 'svg')#one_pair_

    return None
#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
top_metrics = pd.read_csv(args.top_metrics[0])
avdf = pd.read_csv(args.avdf[0])
get_one = bool(args.get_one[0])

cardinality = '_AA20'
#set random seed
np.random.seed(42)
if get_one == True:
    #get one pair per H-group from hgroupdf
    groups = [*Counter(hgroupdf['group']).keys()]
    one_pair_df = pd.DataFrame(columns = hgroupdf.columns)
    for g in groups:
        partial_df = hgroupdf[hgroupdf['group']==g]
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

#Look into outliers
#outliers(catdf, 'pos', outdir)
#outliers(catdf, 'neg', outdir)
#Plot by RCO
# mldists=avdf['ML  distance']
# scores=avdf['lddt_scores_straln']
# plot_rco(catdf, '_straln', 'coolwarm', '1', mldists, scores) #bwr quite good also
# plot_rco(catdf, '_straln', 'coolwarm', '2',mldists, scores) #bwr quite good also

#sig_and_rco(catdf)
#partial_df = catdf[catdf['MLAAdist_straln']>=6]
#partial_df = partial_df[partial_df['MLAAdist_straln']<=8.9]
#plot_class_distr(partial_df, outdir)
#plot_sig(catdf, top_metrics)
#RCO vs deviation
RCO_vs_deviation(catdf, avdf, outdir, '1')
RCO_vs_deviation(catdf, avdf, outdir, '2')
pdb.set_trace()
