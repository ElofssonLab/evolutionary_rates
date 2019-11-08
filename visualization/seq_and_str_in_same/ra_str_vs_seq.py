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

parser.add_argument('--df', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')


def ra_different(df, aln_type, score, cardinality, calc, ylim, outdir):
    '''Produce running average plots for df
    '''

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
    mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
    scores = np.asarray(df[score+aln_type])

    fig = plt.figure(figsize=(10,10)) #set figsize
    for j in np.arange(min(mldists)+step,max(mldists)+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'mean':
            av= np.mean(cut_scores)
        avs.append(av)
        js.append(j-step/2)
        total_avs[j-step] = av
        perc_points.append(len(below_df)/len(df)*100)

        #Include derivatives
        grads = np.gradient(avs)
        gradients[i] = grads
    #Plot RA
    plt.plot(js, avs, linewidth = 1)
    sizes = [js, perc_points] #Save percentage of points in each step
    sns.kdeplot(mldists, scores,  shade=True, shade_lowest = False, cmap = 'Blues')
    plt.title(score)
    plt.xlabel(xlabel)
    plt.ylabel(score)
    plt.ylim(ylim)
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
    fig.savefig(outdir+'running_'+calc+'_'+score+'_'+cardinality+'_'+'.svg', format = 'svg')


    #Plot gradients
    fig = plt.figure(figsize=(10,10)) #set figsize
    plt.scatter(sizes[i][0], gradients[i],s=2)
    plt.ylabel('gradient')
    plt.ylim(grad_ylims[score])
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
    fig.savefig(outdir+'gradient_running_'+calc+'_'+score+'_'+cardinality+'_'+'.svg', format = 'svg')
    #Plot Point distribution
    fig = plt.figure(figsize=(10,10)) #set figsize
    plt.plot(sizes[i][0], sizes[i][1], linewidth = 1)
    plt.xlabel(xlabel)
    plt.ylabel('% of points')
    plt.legend(loc = 'best')
    plt.ylim([0,20])
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
    fig.savefig(outdir+'perc_points_running_'+calc+'_'+score+'_'+cardinality+'_'+'.svg', format = 'svg')

    return None


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
calc = args.calc[0]


aln_types = ['_seqaln', '_straln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.5,1.0]}

cardinality = '_AA20'
for score in ['RMSD','DIFFSS', 'DIFF_ACC', 'lddt_scores']:
    for aln_type in aln_types:
        ylim = ylims[score]
        ra_different(df, aln_type, score, cardinality, calc, ylim, outdir)
