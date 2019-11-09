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


def ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir):
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
    step = 0.1
    top_mldists = np.asarray(topdf['MLAAdist'+cardinality+aln_type])
    top_scores = np.asarray(topdf[score+aln_type])
    hgroup_mldists = np.asarray(hgroupdf['MLAAdist'+cardinality+aln_type])
    hgroup_scores = np.asarray(hgroupdf[score+aln_type])
    
    mldists = np.append(top_mldists, hgroup_mldists)
    scores = np.append(top_scores, hgroup_scores)
    df = pd.concat([topdf, hgroupdf])
    for j in np.arange(min(mldists)+step,max(mldists)+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)
        js.append(j-step/2)
        perc_points.append(len(below_df)/len(df)*100)

    #Include derivatives
    gradients = np.gradient(avs)


    return js, avs, perc_points, mldists, scores, gradients  
   
def make_plots(results, cardinality, outdir, suffix):
    '''Make the plots'''
    
    classes = {1.:'Alpha', 2.: 'Beta', 3.: 'Alpha Beta', 4.: 'Few SS'}
    colors = {1.: 'royalblue', 2.: 'k', 3.: 'green', 4.: 'violet'}
    cmaps = {1.: 'Blues', 2.: 'Greys', 3.: 'Greens', 4.: 'Purples'}
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025]}
    plt.rc('axes', titlesize=10) #set title and label text sizes
    fig = plt.figure(figsize=(10,10)) #set figsize

    for i in [1.,2.,3.,4.]:
        js, avs, perc_points, mldists, scores, gradients = results[i]
        #Plot RA
        plt.plot(js, avs, linewidth = 2, c = colors[i], label = classes[i])
        sns.kdeplot(mldists, scores,  shade=True, shade_lowest = False, cmap = cmaps[i])
        #plt.scatter(mldists, scores, s = 1, c = colors[i], alpha = 0.5, label = classes[i])
        plt.xlabel(xlabel)
        plt.ylabel(score)
        plt.legend()
        plt.ylim(ylim)
        plt.xlim([0,10])
        plt.xticks([0,1,2,3,4,5,6,7,8,9,10])
    fig.savefig(outdir+'class_running_'+suffix, format = 'svg')
    plt.close()

    return None


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]

aln_types = ['_seqaln', '_straln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0]}
cardinality = '_AA20'
results = {}
pdb.set_trace()
for score in ['lddt_scores']:#['RMSD','DIFFSS', 'DIFF_ACC', 'lddt_scores']:
    ylim = ylims[score]
    for aln_type in aln_types:
        for C in [1.,2.,3.,4.]: #Plot per class
            df1 = topdf[topdf['C._x']==C]
            df2 = hgroupdf[hgroupdf['C._x']==C]
            js, avs, perc_points, mldists, scores, gradients = ra_different(df1, df2, aln_type, score, cardinality, calc, ylim, outdir)
            results[C] = [js, avs, perc_points, mldists, scores, gradients]

        suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
        make_plots(results, cardinality, outdir, suffix)
