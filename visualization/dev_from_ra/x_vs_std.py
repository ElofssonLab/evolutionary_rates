#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy.stats import pearsonr
import time

#Custom
from percent_ss import parse_ss
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

def AA6_distribution(df, aln_type):
    '''Calculate AA6 distribution in sequence (gaps have their own state as well)
    Groupings: K=[KR],D=[DE],Y=[YWFH],T=[TSQN],C=[CVMLIA], P=[PG], "-"="-"
    '''
    #Save percentages:
    percentage_dict = {
    'K1'+aln_type:[],
    'K2'+aln_type:[],
    'D1'+aln_type:[],
    'D2'+aln_type:[],
    'Y1'+aln_type:[],
    'Y2'+aln_type:[],
    'T1'+aln_type:[],
    'T2'+aln_type:[],
    'C1'+aln_type:[],
    'C2'+aln_type:[],
    'P1'+aln_type:[],
    'P2'+aln_type:[],
    '-1'+aln_type:[],
    '-2'+aln_type:[]
    }


    #reset index
    df = df.reset_index()
    t1 = time.time()
    for i in range(len(df)):
        row = df.iloc[i]
        l1 = len(row['AA6_1_straln'])
        l2 = len(row['AA6_2_straln'])
        s1 = Counter(row['AA6_1_straln'])
        s2 = Counter(row['AA6_2_straln'])

        #Get s1 ad s2 percentages
        for key in ['P', 'C', 'K', 'T', 'D', 'Y', '-']:
            try:
                percentage_dict[key+'1'+aln_type].append(s1[key]/l1)
            except:
                percentage_dict[key+'1'+aln_type].append(0)
            try:
                percentage_dict[key+'2'+aln_type].append(s2[key]/l2)
            except:
                percentage_dict[key+'2'+aln_type].append(0)

    print('Calculating AA6',time.time()-t1,'s')

    for key in percentage_dict:
        df[key] = percentage_dict[key]

    return df, percentage_dict.keys()

def dev_from_av(avdf, df, score, aln_type, cardinality):
    '''Calculate number of stds from total running average for each pair
    '''
    start = 0
    end = 6
    step = 0.1
    std_df = pd.DataFrame() #save vals
    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_df = df[df['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]

        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])

        #Get total average in interval
        x = np.round(j-step/2, 2)
        tav = avdf[avdf['ML  distance']==x][score+aln_type].values[0] #total average in current interval

        #Deviation
        dev = cut_scores-tav
        below_df[score+aln_type+'_av'] = [tav]*len(dev)
        below_df[score+aln_type+'_dev'] = dev
        std_df = pd.concat([std_df, below_df])

    all_devs = np.absolute(np.array(std_df[score+aln_type+'_dev']))
    n_within = np.where(all_devs<0.05)[0].size
    print('Fraction of points within 0.05 lddt fromtotal ra:', str(n_within)+'/'+str(len(std_df)),n_within/len(std_df))
    return std_df

def plot_format(fig, ax, outname, ylabel):
    '''Format plot layout
    '''

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(outname, format = 'png')
    plt.close()

def plot_x_vs_std(std_df, avdf, single_features, double_features, score, aln_type, outdir):
    '''Plot different features agains the std deviation for each pair
    '''
    #Make outdir
    try:
        os.mkdir(outdir)
    except:
        print('Directory '+outdir+' exists')
    matplotlib.rcParams.update({'font.size': 7})
    #Plot all pairs
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54))
    ax.scatter(std_df['MLAAdist'+aln_type], std_df[score+aln_type], s= 0.1, c='cornflowerblue', label = 'Pair' ) #All points
    sns.kdeplot(std_df['MLAAdist'+aln_type], std_df[score+aln_type], shade = False, cmap = 'Blues')
    ax.plot(avdf['ML  distance'], avdf[score+aln_type], color = 'darkblue', linewidth = 1, label = 'Running average')
    #plot intercept
    ax.plot(avdf['ML  distance'], np.array(avdf[score+aln_type])+0.05, '--', c = 'darkblue', linewidth = 1) #positive stds
    ax.plot(avdf['ML  distance'], np.array(avdf[score+aln_type])-0.05, '--', c = 'darkblue', linewidth = 1, label = '+/- (1-intercept)') #negative stds
    ax.set_xlim([0,6.1])
    ax.set_xticks([0,1,2,3,4,5,6])
    ax.set_xlabel('AA20 ED')
    ax.legend(frameon=False, markerscale = 10)
    if score == 'lddt_scores':
        ylabel = 'lDDT score'
    else:
        ylabel = score
    plot_format(fig, ax, outdir+'all_pairs.png', ylabel)

    #Plot sequence and structure
    if score+aln_type == 'lddt_scores_straln':
        fig, ax = plt.subplots(figsize=(9/2.54,9/2.54))
        ax.plot(avdf['ML  distance'], avdf[score+aln_type], color = 'darkblue', linewidth = 1, label = 'Structure RA')
        sns.kdeplot(std_df['MLAAdist'+aln_type], std_df[score+aln_type], shade = False, cmap = 'Blues')
        ax.plot(avdf['ML  distance'], avdf[score+'_seqaln'], color = 'darkgreen', linewidth = 1, label = 'Sequence RA')
        sns.kdeplot(std_df['MLAAdist_seqaln'], std_df[score+'_seqaln'], shade = False, cmap = 'Greens')
        ax.set_xlim([0,6.1])
        ax.set_xticks([0,1,2,3,4,5,6])
        ax.set_ylim([0.2,1])
        ax.set_xlabel('AA20 ED')
        ax.legend(frameon=False)
        plot_format(fig, ax, outdir+'seq_vs_str.png', ylabel)
        print(score,pearsonr(std_df[score+'_seqaln'], std_df[score+aln_type]))
        print('ED',pearsonr(std_df['MLAAdist_seqaln'], std_df['MLAAdist_straln']))

    #Plot hist of distances
    fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
    sns.distplot(std_df[score+aln_type+'_dev'], hist = True)
    ax.set_xlabel('Deviation from line')
    plot_format(fig, ax, outdir+'all_pairs_hist.png', 'Density')

    #Save pearsonr
    p_R = {}
    titles = {'RCO':'RCO', 'aln_len'+aln_type:'Aligned length', 'l':'Length', 'percent_aligned'+aln_type:'% Aligned',
    'K':'KR','D':'DE','Y':'YWFH','T':'TSQN','C':'CVMLIA', 'P':'PG', '-':'Gap', 'L':'Loop', 'S':'Sheet', 'H':'Helix',
    score+aln_type+'classes':'Class', score+aln_type+'_sizes':'Group size', 'CD': 'Contact Density', 'MLAAdist'+aln_type: 'AA20 ED'}

    #Color in outliers
    pos_odf = std_df[std_df[score+aln_type]>0.9]
    pos_odf = pos_odf[pos_odf['MLAAdist'+aln_type]>2]
    counts = Counter(pos_odf['group'])

    neg_odf = std_df[std_df[score+aln_type]<0.45]
    neg_odf = neg_odf[neg_odf['MLAAdist'+aln_type]<2]
    counts = Counter(neg_odf['group'])

    for x in single_features:
        fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
        plt.scatter(std_df[x], std_df[score+aln_type+'_dev'], label = x,s=0.3, color = '#1f77b4', alpha = 0.2)
        sns.kdeplot(std_df[x], std_df[score+aln_type+'_dev'], shade = False, cmap = 'Blues')
        #Plot outlier group
        plt.scatter(pos_odf[x], pos_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'maroon')
        plt.scatter(neg_odf[x], neg_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'k')
        ax.set_xlabel(titles[x])
        plot_format(fig, ax, outdir+x+'.png', 'Deviation from line')
        p_R[titles[x]] = pearsonr(std_df[x], std_df[score+aln_type+'_dev'])[0] #returns (Pearson’s correlation coefficient, 2-tailed p-value)

    #Double
    for x in double_features:
        fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
        if x == 'RCO' or x == 'CD':
            for i in ['1', '2']:
                plt.scatter(std_df[x+i], std_df[score+aln_type+'_dev'], s=0.3, color = '#1f77b4', label = 'Domain 1', alpha = 0.2)
                sns.kdeplot(std_df[x+i], std_df[score+aln_type+'_dev'], shade = False, cmap = 'Blues')
                plt.scatter(pos_odf[x+i], pos_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'maroon')
                plt.scatter(neg_odf[x+i], neg_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'k')
                ax.set_xlabel(titles[x]+i)
                plot_format(fig, ax, outdir+x+i+'.png', 'Deviation from line')
                p_R[titles[x]+i] = pearsonr(std_df[x+i], std_df[score+aln_type+'_dev'])[0] #returns (Pearson’s correlation coefficient, 2-tailed p-value)

        else:
            for i in ['1', '2']:
                plt.scatter(std_df[x+i+aln_type], std_df[score+aln_type+'_dev'] ,s=0.3, color = '#1f77b4', label = 'Domain 1', alpha = 0.2)
                sns.kdeplot(std_df[x+i+aln_type], std_df[score+aln_type+'_dev'], shade = False, cmap = 'Blues')
                #Plot outlier group
                plt.scatter(pos_odf[x+i+aln_type], pos_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'maroon')
                plt.scatter(neg_odf[x+i+aln_type], neg_odf[score+aln_type+'_dev'], label = x,s=0.3, color = 'k')
                ax.set_xlabel(titles[x]+i)
                plot_format(fig, ax, outdir+titles[x]+i+'.png', 'Deviation from line')
                p_R[titles[x]+i] = pearsonr(std_df[x+i+aln_type], std_df[score+aln_type+'_dev'])[0] #returns (Pearson’s correlation coefficient, 2-tailed p-value)


    #Plot personr
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54))
    pearsondf = pd.DataFrame()
    pearsondf['Feature'] = [*p_R.keys()]
    pearsondf['Pearson R'] = [*p_R.values()]
    pearsondf=pearsondf.sort_values(by='Pearson R',ascending=True)
    sns.barplot(x="Pearson R", y="Feature", data=pearsondf)
    ax.set_xticks([-0.2,0,0.2])
    plot_format(fig, ax, outdir+'pearsonr.png', 'Feature')

    return None

def get_over_10(df):
    '''Get all topologies with at least 10 entries in the interval
    '''

    #Get topologies with 10 points
    topcounts = Counter(df['group'])
    vals = np.array([*topcounts.values()])
    num_tops_with10 = len(np.where(vals>9)[0]) #Get all topologies with at least 10 values
    print('Fraction of topologies with at least 10 entries: '+str(num_tops_with10)+'/'+str(len(vals)),num_tops_with10/len(vals))
    topologies = np.array([*topcounts.keys()])
    topologies = topologies[np.where(vals>9)[0]]

    df10 = pd.DataFrame()
    for top in topologies:
        sel = df[df['group']==top]
        df10 = pd.concat([df10, sel])

    print('Fraction of points with 10: '+str(len(df10))+'/'+str(len(df)),len(df10)/len(df))
    return df10

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

#rename TMscore cols
hgroupdf = hgroupdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})
topdf = topdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})

catdf = pd.concat([topdf, hgroupdf])
#Rename class column
catdf = catdf.rename(columns={'C._x':'Class'})

#The ones should actually be zeros
catdf['RCO1']=catdf['RCO1'].replace([1], 0)
catdf['RCO2']=catdf['RCO2'].replace([1], 0)


cardinality = '' #AA20

for score in ['lddt_scores']:
    for aln_type in ['_straln']:
        #select below 6 using seq or str
        catdf_s = catdf[catdf['MLAAdist'+aln_type]<=6]
        print('Fraction of points under ED 6: '+str(len(catdf_s))+'/'+str(len(catdf)),len(catdf_s)/len(catdf))
        #catdf_s = get_over_10(catdf_s)

        single_features = ['percent_aligned'+aln_type, 'MLAAdist'+aln_type]
        double_features = ['CD', 'RCO', 'l', 'P', 'C', 'K', 'T', 'D', 'Y', '-', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix, contact density
        catdf_s, perc_keys = AA6_distribution(catdf_s, aln_type) #Get AA6 frequencies
        catdf_s = parse_ss(catdf_s, aln_type) #Get % ss
        std_df = dev_from_av(avdf, catdf_s, score, aln_type, cardinality)
        std_df.to_csv('top10df.csv')
        plot_x_vs_std(std_df, avdf, single_features, double_features, score, aln_type, outdir+score+aln_type+'/')
