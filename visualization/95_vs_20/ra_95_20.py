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

parser.add_argument('--df95', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with the 95 % reduction.')

parser.add_argument('--df20', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with the 20 % reduction.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either mean or average.')

#FUNCTIONS
def ra_different(dfs, aln_type, score, cardinality, calc, plot_num, pdf, fig, ylim, title):
    '''Produce running average plots for df
    '''

    colors = {0: 'b', 1: 'darkred'}
    labels = {0:'95% reduction', 1: '20% reduction'}
    titles = {'seqaln': 'Sequence alignments', 'straln': 'Structure alignments'}
    plt.rc('axes', titlesize=10) #set title and label text sizes
    plt.subplot(plot_num) #set plot_num
    sizes = {} #Save percentage of points in each step
    xlabel = 'ML '+cardinality[1:]+' distance'
    title = titles[aln_type[1:]]
    gradients = {}
    for i in range(2):
        df = dfs[i]
        #Plot total average for cardinality
        color = colors[i]
        if cardinality == '_AA20':
            cardinality = ''
        avs = [] #Save average score
        js = [] #Save dists
        perc_points = []
        total_avs = {}
        step = 0.1
        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score+aln_type])

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
        plt.plot(js, avs, label = labels[i], linewidth = 1, color = colors[i])
        sizes[i] = [js, perc_points]

        #plt.scatter(mldists, scores, color = color, s= 1)
    plt.legend(loc = 'best')
    #plt.xlabel('MLAAdist'+aln_type)
    plt.ylabel(score)
    plt.ylim(ylim)
    plt.xlim([0,6])
    plt.title(title)

    #Plot gradients
    plot_num+=3
    plt.subplot(plot_num) #set plot_num
    for i in range(2):
        plt.scatter(sizes[i][0], gradients[i],s=2, label = labels[i], color = colors[i])
    plt.ylabel('gradient')
    plt.ylim([-0.1,0.1])
    plt.xlim([0,6])
    plt.legend(loc = 'best')
    #Plot Point distribution
    plot_num += 3
    plt.subplot(plot_num) #set plot_num
    for i in range(2):
        plt.plot(sizes[i][0], sizes[i][1], label = labels[i], linewidth = 1,  color = colors[i])
    plt.xlabel(xlabel)
    plt.ylabel('% of points')
    plt.legend(loc = 'best')
    plt.ylim([0,20])
    plt.xlim([0,6])

    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df95 = pd.read_csv(args.df95[0])
df20 = pd.read_csv(args.df20[0])
outdir = args.outdir[0]
calc = args.calc[0]
score ='RMSD'
cardinality = '_AA20'
pdf = PdfPages(outdir+score+'_'+calc+'_95_vs_20.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize

aln_types = ['_seqaln', '_straln']
title = 'Sequence vs Structure alignments'
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.5, 1.0]}



plot_num = 331
ylim = ylims[score]
dfs = [df95, df20]
for aln_type in aln_types:
    pdf, fig = ra_different(dfs, aln_type, score, cardinality, calc, plot_num, pdf, fig, ylim, title)
    plot_num +=2

pdf.savefig(fig)
fig.savefig(outdir+score+'_'+calc+cardinality+'_95_vs_20.svg', format = 'svg')
pdf.close()
