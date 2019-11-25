#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy import stats
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that calculates deviations for running averages.''')

parser.add_argument('--topdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')

parser.add_argument('--avdf', nargs=1, type= str,
default=sys.stdin, help = 'Dataframe with averages.')

def dev_from_av(avdf, df, score, aln_type, cardinality, max_seqdist):
    '''Calculate avearge dev from total running average within group and significance
    '''

    if cardinality == '_AA20':
            cardinality = ''

    absolute_deviations = [] #absolute of all devs
    rco_extreme = []
    avs = [] #Save average score
    avdevs = [] #Save deviations
    step = 0.1


    mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
    start = np.round(min(mldists),2)
    start = float(str(start)[0:3])
    end = min(max(mldists), max_seqdist) #End at max_seqdist
    scores = np.asarray(df[score+aln_type])
    js = [] #Save js

    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_df = df[df['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]

        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        if len(below_df)<1: #May be discontinuous in step
            continue
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)

        x = np.round(j-step/2, 2)
        js.append(x)
        tav = avdf[avdf['ML  distance']==x][score+aln_type].values[0] #total average in current interval
        avdevs.append(av-tav)

        #Plot deviation against RCO
        absolute_deviations.extend(np.absolute(cut_scores-tav))
        rco1 = [*below_df['RCO1']]
        rco2 = [*below_df['RCO2']]


    #Do a t-test to calculate if the deviation is significant for each H-group
    #The variances will be unequal, as each H-group only encompasses very feq points
    #True values = 0, no deviation from total average
    truevals = np.zeros(len(df)) #size should not matter since zeros
    statistic, pvalue = stats.ttest_ind(avdevs, truevals, equal_var = False)
    #plt.plot(js, avs)

    return np.average(avdevs), pvalue, js, avs

def plot_partial(partial_df, avdf, catdf, name):
    '''RA plots of partial dfs and the total RA
    '''
    topologies = [*partial_df['Topology']]
    mldists = [*partial_df['lddt_scores_straln_seqdists']]
    scores = [*partial_df['lddt_scores_straln_ra']]
    fig = plt.figure(figsize=(10,10)) #set figsize
    matplotlib.rcParams.update({'font.size': 22})
    for i in range(len(partial_df)):
        top = topologies[i]
        plt.plot(mldists[i],scores[i], alpha = 0.1, color = 'b', linewidth =2)
    plt.plot(avdf['ML  distance'], avdf['lddt_scores_straln'], color = 'r', linewidth = 3, label = 'Total RA')
    plt.legend()
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.ylim([0.2,1])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score for straln workflow')
    fig.savefig(outdir+name, format = 'png')

    #Scatterplot
    fig = plt.figure(figsize=(10,10)) #set figsize
    for i in range(len(partial_df)):
        top = topologies[i]
        topdf = catdf[catdf['group']==top]
        plt.scatter(topdf['MLAAdist_straln'],topdf['lddt_scores_straln'],alpha = 0.2, color = 'b', s = 1)
    plt.plot(avdf['ML  distance'], avdf['lddt_scores_straln'], color = 'r', linewidth = 3, label = 'Total RA')
    plt.legend()
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    plt.ylim([0.2,1])
    plt.xlabel('ML AA20 distance')
    plt.ylabel('lddt score for straln workflow')
    fig.savefig(outdir+'scatter_'+name, format = 'png')

    return None

def class_percentages(df):
    '''Calculate class percentages
    '''

    percentages = []
    for C in [1.,2.,3.,4.]:
        try:
            percentages.append(len(df[df['C._x']==C])/len(df))
        except:
            percentages.append(0)

    return percentages

def anova(cat_dev):
    '''Perform anova
    '''

    #Make violinplots
    features = ['RCO1', 'RCO2', 'aln_len_straln', 'l1_straln', 'l2_straln', 'percent_aligned_straln']
    for feature in features:
        matplotlib.rcParams.update({'font.size': 22})
        fig = plt.figure(figsize=(10,10)) #set figsize
        sns.violinplot(data = cat_dev, x = 'Significance', y = feature)
        fig.savefig(outdir+feature+aln_type+score+'.png', format = 'png')
    #Calculate fraction retained
    eta_sq = []
    omega_sq = []
    Pr = []

    #Fit ANOVA
    for column in features:
        test_str = column+' ~ C(Significance)'
        results = ols(test_str, data=cat_dev).fit()
        aov_table = sm.stats.anova_lm(results, typ=2)
        outp = anova_table(aov_table)
        eta_sq.append(outp.iloc[0]['eta_sq'])
        omega_sq.append(outp.iloc[0]['omega_sq'])
        Pr.append(outp.iloc[0]['PR(>F)'])

    anova_df =  pd.DataFrame(data = {'Feature': features, 'eta_sq': eta_sq, 'omega_sq':omega_sq, 'PR(>F)':Pr})

    return anova_df

def anova_table(aov):
    aov['mean_sq'] = aov[:]['sum_sq']/aov[:]['df']

    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])

    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*aov['mean_sq'][-1]))/(sum(aov['sum_sq'])+aov['mean_sq'][-1])

    cols = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[cols]
    return aov

def ttest_table(neg_sig, pos_sig, nonsig_df, features):
    '''Perform t-tests for features
    '''

    print('Feature\tNeg\tPositive')
    for feature in features:
        print(feature+'\t')
        statistic, pvalue = stats.ttest_ind(neg_sig[feature], nonsig_df[feature], equal_var = False)
        print(str(statistic)+'\t')
        statistic, pvalue = stats.ttest_ind(pos_sig[feature], nonsig_df[feature], equal_var = False)
        print(str(statistic)+'\n')


    return None
#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
avdf = pd.read_csv(args.avdf[0])

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
catdf = pd.concat([topdf, hgroupdf])
#select below 6
catdf = catdf[catdf['MLAAdist_straln']<=6]
topcounts = Counter(catdf['group'])
vals = np.array([*topcounts.values()])
num_tops_with10 = len(np.where(vals>9)[0]) #Get all topologies with at least 10 values
print('Fraction of topologies with at least 10 entries:',num_tops_with10/len(vals))
topologies = np.array([*topcounts.keys()])
topologies = topologies[np.where(vals>9)[0]]


#Save pvalues
top_metrics = pd.DataFrame()
top_metrics['Topology'] = topologies

for score in ['lddt_scores']:
    for aln_type in ['_straln']:
        avs_from_line = [] #save avs from line and pvals
        pvals = []
        all_js = []
        all_avs = []
        for top in topologies:
            df = catdf[catdf['group']==top]
            av_from_line, pvalue, js, avs = dev_from_av(avdf, df, score, aln_type, cardinality, 6)
            avs_from_line.append(av_from_line)
            pvals.append(pvalue)
            all_js.append(js)
            all_avs.append(avs)
        top_metrics[score+aln_type+'_pval'] = pvals
        top_metrics[score+aln_type+'_av_dev'] = avs_from_line
        top_metrics[score+aln_type+'_seqdists'] = all_js
        top_metrics[score+aln_type+'_ra'] = all_avs


#Get significant
sig_df = top_metrics[top_metrics['lddt_scores_straln_pval']<0.05/len(top_metrics)]
#Get pos deviating>0
pos_sig = sig_df[sig_df['lddt_scores_straln_av_dev']>0]
pos_sig['Significance'] = 'Positive'
#Get neg deviating<0
neg_sig = sig_df[sig_df['lddt_scores_straln_av_dev']<0]
neg_sig['Significance'] = 'Negative'
#check
if len(pos_sig)+len(neg_sig) != len(sig_df):
    pdb.set_trace()

#Plot the RAs of the pos and neg sig groups
#plot_partial(pos_sig, avdf, catdf, 'ra_pos_sig.png')
#plot_partial(neg_sig, avdf, catdf, 'ra_neg_sig.png')


#Get data from cat_df matching sig topologies
pos_sig = pd.merge(pos_sig, catdf, left_on='Topology', right_on='group', how='left')
neg_sig = pd.merge(neg_sig, catdf, left_on='Topology', right_on='group', how='left')
#Get non significant
nonsig_df = top_metrics[top_metrics['lddt_scores_straln_pval']>=0.05/len(top_metrics)]
#plot_partial(nonsig_df, avdf, catdf, 'ra_non_sig.png')
#Get data from cat_df matching sig topologies
nonsig_df = pd.merge(nonsig_df, catdf, left_on='Topology', right_on='group', how='left')
nonsig_df['Significance'] = 'Non-significant'

#Concat
cat_dev = pd.concat([pos_sig, neg_sig])
classes = ['Mainly Alpha', 'Mainly Beta', 'Alpha Beta', 'Few SS']
print('%pos sig', class_percentages(pos_sig))
print('%neg sig', class_percentages(neg_sig))
print('%nonsig', class_percentages(nonsig_df))
cat_dev = pd.concat([cat_dev, nonsig_df])
#Fraction of pairs retained
print('Fraction of pairs within topologies with at least 10 entries:', len(cat_dev)/len(catdf))

#Perform t-tests
features = ['RCO1', 'RCO2', 'aln_len_straln', 'l1_straln', 'l2_straln']
ttest_table(neg_sig, pos_sig, nonsig_df, features)
#Calculate ANOVA
#anova(cat_dev)

#plot_partial(sig_df)
top_metrics.to_csv(outdir+'top_metrics.csv')
pdb.set_trace()
