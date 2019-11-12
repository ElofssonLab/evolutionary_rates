#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots running averages.''')

parser.add_argument('--topdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')

parser.add_argument('--get_one', nargs=1, type= int,
default=sys.stdin, help = 'Get one pair from each H-group (1) or all (0).')


def ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir, av_df):
    '''Produce running average plots for df
    '''

    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    colors = {'_seqaln': 'k', '_straln': 'r'}
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025]}
    plt.rc('axes', titlesize=10) #set title and label text sizes

    gradients = {}

    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''
    avs = [] #Save average score
    js = [] #Save dists
    perc_points = []
    total_avs = {}
    step = 0.1
    top_mldists = np.asarray(topdf['MLAAdist'+cardinality+aln_type])
    top_scores = np.asarray(topdf[score+aln_type])
    hgroup_mldists = np.asarray(hgroupdf['MLAAdist'+cardinality+aln_type])
    hgroup_scores = np.asarray(hgroupdf[score+aln_type])

    mldists = np.append(top_mldists, hgroup_mldists)
    scores = np.append(top_scores, hgroup_scores)
    df = pd.concat([topdf, hgroupdf])
    fig = plt.figure(figsize=(10,10)) #set figsize
    for j in np.arange(min(mldists)+step,max(mldists)+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)
        js.append(np.round(j-step/2,2))
        total_avs[j-step] = av
        perc_points.append(len(below_df)/len(df)*100)

    #Include derivatives
    gradients = np.gradient(avs)
    #Plot RA
    plt.plot(js, avs, linewidth = 2, c = 'b', label = 'Running average')
    #sns.kdeplot(mldists, scores,  shade=True, shade_lowest = False, cmap = 'Blues')
    plt.scatter(top_mldists, top_scores, s = 1, c = 'k', alpha = 1.0, label = 'topology')
    plt.scatter(hgroup_mldists, hgroup_scores, s = 1, c = 'r', alpha = 1.0, label = 'hgroups')
    plt.xlabel(xlabel)
    plt.ylabel(score)
    plt.legend()
    plt.ylim(ylim)
    plt.xlim([0,10])
    plt.xticks([0,1,2,3,4,5,6,7,8,9,10])
    fig.savefig(outdir+'running_'+suffix, format = 'svg')
    plt.close()

    #Plot gradients
    fig = plt.figure(figsize=(10,10)) #set figsize
    plt.scatter(js, gradients,s=5)
    plt.plot(js, gradients, linewidth = 1)
    plt.ylabel('gradient')
    #plt.ylim(grad_ylims[score])
    #plt.xlim([2,10])
    #plt.xticks([2,3,4,5,6,7,8,9,10])
    fig.savefig(outdir+'gradient_running_'+suffix, format = 'svg')
    #Plot Point distribution
    fig = plt.figure(figsize=(10,10)) #set figsize
    plt.plot(js, perc_points, linewidth = 2)
    plt.xlabel(xlabel)
    plt.ylabel('% of points')
    #plt.xlim([2,10])
    #plt.xticks([2,3,4,5,6,7,8,9,10])
    fig.savefig(outdir+'perc_points_running_'+suffix, format = 'svg')

    av_df['ML '+cardinality[1:]+' distance'] = js
    av_df[score+aln_type] = avs
    return av_df


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
get_one = bool(args.get_one[0])

aln_types = ['_seqaln', '_straln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0]}

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

cardinality = '_AA20'
av_df = pd.DataFrame()
for score in ['RMSD','DIFFSS', 'DIFF_ACC', 'lddt_scores']:
    for aln_type in aln_types:
        ylim = ylims[score]
        av_df = ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir, av_df)
av_df.to_csv(outdir+'av_df.csv')
