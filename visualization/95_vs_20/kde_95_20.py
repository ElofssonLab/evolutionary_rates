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
parser = argparse.ArgumentParser(description = '''A program that plots KDE plots.''')

parser.add_argument('--df95', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with the 95 % reduction.')

parser.add_argument('--df20', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with the 20 % reduction.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either mean or average.')

#FUNCTIONS
def kde(dfs, aln_type, score, cardinality, calc, plot_num, pdf, fig, ylim):
    '''Produce running kde plots for df
    '''

    #Start with kde plot
    cmaps = {0: 'Blues', 1: 'Oranges'}
    colors = {0: 'b', 1: 'r'}
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

        color = colors[i]
        if cardinality == '_AA20':
            cardinality = ''

        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score+aln_type])

        #Plot KDE
        sns.kdeplot(mldists, scores, label = labels[i], cmap = cmaps[i], shade = True, shade_lowest=False)



    plt.legend(loc = 'best')
    plt.ylabel(score)
    plt.ylim(ylim)
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
    plt.title(title)


    #Do individual kde plots as well
    for i in range(2):
        df = dfs[i]
        plot_num +=3
        plt.subplot(plot_num) #set plot_num

        color = colors[i]
        if cardinality == '_AA20':
            cardinality = ''


        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score+aln_type])
        sns.kdeplot(mldists, scores, cmap = cmaps[i], shade = True, shade_lowest=False)
        plt.ylabel(score)
        plt.ylim(ylim)
        plt.xlim([0,6])
        plt.xticks([0,1,2,3,4,5,6])
        plt.title(labels[i])
        if i == 1:
            plt.xlabel(xlabel)


    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df95 = pd.read_csv(args.df95[0])
df20 = pd.read_csv(args.df20[0])
outdir = args.outdir[0]
calc = args.calc[0]
score ='RMSD'
cardinality = '_AA20'
pdf = PdfPages(outdir+score+'_kde_95_vs_20.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize

aln_types = ['_seqaln', '_straln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6]}



plot_num = 331
ylim = ylims[score]
dfs = [df95, df20]
for aln_type in aln_types:
    pdf, fig = kde(dfs, aln_type, score, cardinality, calc, plot_num, pdf, fig, ylim)
    plot_num +=2

pdf.savefig(fig)
fig.savefig(outdir+score+'_kde_95_vs_20.svg', format = 'svg')
pdf.close()
