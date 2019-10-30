#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
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
def ra_different(df, aln_type, score, cardinalities, calc, plot_num, pdf, fig, ylim, title):
    '''Produce running average plots for df
    '''

    viridis = cm.get_cmap('viridis', 12)
    colors = {'_AA2':viridis(0.8), '_AA3':viridis(0.6), '_AA6':viridis(0.45), '_AA20':viridis(0.1)} #Set colors

    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025]}
    plt.rc('axes', titlesize=10) #set title and label text sizes
    plt.subplot(plot_num) #set plot_num
    sizes = {} #Save percentage of points in each step

    gradients = {}

    for cardinality in cardinalities:
        label = cardinality[1:] #set label
        color = colors[cardinality]
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
            if calc == 'median':
                av= np.median(cut_scores)
            avs.append(av)
            js.append(j-step/2)
            total_avs[j-step] = av
            perc_points.append(len(below_df)/len(df)*100)

        #Include derivatives
        grads = np.gradient(avs)
        gradients[cardinality] = grads
        #Plot RA
        plt.plot(js, avs, label = label, linewidth = 1, color = color)
        sizes[cardinality] = [js, perc_points] #save the mlaadists and perc points

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
    for cardinality in cardinalities:
        label = cardinality[1:] #set label
        color = colors[cardinality] #set color
        if cardinality == '_AA20':
            cardinality = ''
        plt.scatter(sizes[cardinality][0], gradients[cardinality],s=2, label = label, color = color)
    plt.ylabel('gradient')
    #plt.ylim(grad_ylims[score])
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])
    plt.legend(loc = 'best')
    #Plot Point distribution
    plot_num += 3
    plt.subplot(plot_num) #set plot_num
    for cardinality in cardinalities:
        label = cardinality[1:] #set label
        color = colors[cardinality] #set color
        if cardinality == '_AA20':
            cardinality = ''
        plt.plot(sizes[cardinality][0], sizes[cardinality][1], label = label, linewidth = 1, color = color)
    plt.xlabel('ML distance')
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

cardinalities = ['_AA2', '_AA3', '_AA6', '_AA20']
score ='RMSD'
pdf = PdfPages(outdir+score+'_'+calc+'_all_cards.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize

aln_types = ['_seqaln', '_straln']
title = 'Sequence vs Structure alignments'
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.5,1.0]}



plot_num = 331
for aln_type in aln_types:
    ylim = ylims[score]
    pdf, fig = ra_different(df, aln_type, score, cardinalities, calc, plot_num, pdf, fig, ylim, title)
    plot_num +=2

pdf.savefig(fig)
fig.savefig(outdir+score+'_'+calc+'_all_cards.svg', format = 'svg')
pdf.close()
