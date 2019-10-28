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
default=sys.stdin, help = 'path to directory with pdb files.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either mean or average.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')

def count_groups():
    '''Get the original size for each H-group
    '''
    H_group_sizes = {}
    for uid in H_groups:
        group = H_groups[uid]
        if group not in H_group_sizes.keys():
            H_group_sizes[group] = 1
        else:
            H_group_sizes[group] += 1
def ra_different(df, aln_types, score, cardinality, calc, plot_num, pdf, fig, ylim, title):
    '''Produce running average plots for df
    '''

    colors = {0: 'k', 1: 'r'}
    labels = {0:'Sequence alignments', 1: 'Structure alignments'}
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025]}
    plt.rc('axes', titlesize=10) #set title and label text sizes
    plt.subplot(plot_num) #set plot_num
    sizes = {} #Save percentage of points in each step

    gradients = {}
    for i in range(2):
        aln_type = aln_types[i]
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

    plt.title(score)
    plt.ylabel(score)
    plt.ylim(ylim)
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])



    #Plot gradients
    plot_num+=3
    plt.subplot(plot_num) #set plot_num
    for i in range(2):
        plt.scatter(sizes[i][0], gradients[i],s=2, label = labels[i], color = colors[i])
    plt.ylabel('gradient')
    plt.ylim(grad_ylims[score])
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
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
    plt.xticks([0,1,2,3,4,5,6])

    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
calc = args.calc[0]
plot_gradients = args.plot_gradients[0]

cardinality = '_AA20'
score ='RMSD'
pdf = PdfPages(outdir+'rmsd_lddt_'+calc+'_str_vs_seq.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize

aln_types = ['_seqaln', '_straln']
title = 'Sequence vs Structure alignments'
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.5,1.0]}



plot_num = 331
for score in ['RMSD', 'lddt_scores']:
    ylim = ylims[score]
    pdf, fig = ra_different(df, aln_types, score, cardinality, calc, plot_num, pdf, fig, ylim, title)
    plot_num +=2

pdf.savefig(fig)
fig.savefig(outdir+'rmsd_lddt_'+calc+'_str_vs_seq.svg', format = 'svg')
pdf.close()
