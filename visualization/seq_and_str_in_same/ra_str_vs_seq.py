#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib
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
    '''Produce running average plots for df using 100 steps with an equal amount of points in each
    '''

    #Double column = 183 mm, maxfont = 8
    matplotlib.rcParams.update({'font.size': 7})
    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025]}


    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''
    avs = [] #Save average score
    stds = [] #Save std deviation
    js = [] #Save dists
    perc_points = [] #Save percentage of points within seqdist interval
    perc_within_acc = []
    perc_within_std = []

    top_mldists = np.asarray(topdf['MLAAdist'+cardinality+aln_type])
    top_scores = np.asarray(topdf[score+aln_type])
    hgroup_mldists = np.asarray(hgroupdf['MLAAdist'+cardinality+aln_type])
    hgroup_scores = np.asarray(hgroupdf[score+aln_type])

    mldists = np.append(top_mldists, hgroup_mldists)
    scores = np.append(top_scores, hgroup_scores)
    df = pd.concat([topdf, hgroupdf])
    fig = plt.figure(figsize=(9/2.54,9/2.54)) #set figsize


    #Sort df by x-value
    df = df.sort_values(by=['MLAAdist'+cardinality+aln_type], ascending=True)
    step = 0.1

    mldists = np.append(top_mldists, hgroup_mldists)
    scores = np.append(top_scores, hgroup_scores)

    for j in np.arange(min(mldists)+step,max(mldists)+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)
        std = np.std(cut_scores)
        stds.append(std)
        js.append(np.round(j-step/2,2))
        perc_points.append(len(below_df)/len(df)*100)
        #Get % points within 0.05
        outside = np.where(cut_scores<av-0.05)[0].size+np.where(cut_scores>av+0.05)[0].size
        perc_within_acc.append((cut_scores.size-outside)/cut_scores.size)
        #Get % points within std
        outside_std = np.where(cut_scores<av-std)[0].size+np.where(cut_scores>av+std)[0].size
        perc_within_std.append((cut_scores.size-outside_std)/cut_scores.size)
    #Include derivatives
    gradients = np.gradient(avs)
    #Fit polyline
    z = np.polyfit(js, avs, deg = 3)
    p = np.poly1d(z)
    #Plot RA
    plt.plot(js, avs, linewidth = 1, c = 'g', label = 'Running average')
    #sns.kdeplot(mldists, scores,  shade=True, shade_lowest = False, cmap = 'Blues')
    plt.scatter(hgroup_mldists, hgroup_scores, s = 0.1, c = 'lightseagreen', alpha = 0.5, label = 'Dataset 1')
    plt.scatter(top_mldists, top_scores, s = 0.1, c = 'b', alpha = 1.0, label = 'Dataset 3')
    #plot stddev
    plt.plot(js, np.array(avs)+np.array(stds), '--', c = 'g', linewidth = 1) #positive stds
    plt.plot(js, np.array(avs)-np.array(stds), '--', c = 'g', linewidth = 1, label = 'Standard deviation') #negative stds
    #Plot polynomial fit
    plt.plot(js,p(js), label = '3 dg polynomial fit',linewidth = 1, c= 'indigo')
    plt.xlabel(xlabel)
    if score == 'lddt_scores':
        plt.ylabel('lDDT score')
    else:
        plt.ylabel(score)
    plt.legend(markerscale=5)
    plt.ylim(ylim)
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    #pdb.set_trace()
    fig.savefig(outdir+'running_'+suffix, format = 'svg')
    plt.close()

    fig = plt.figure(figsize=(9/2.54,4.5/2.54))
    #Plot Percent within accurac
    plt.plot(js,np.round(np.array(perc_wihtin_std)*100,2), linewidth = 1, c= 'g',  label = 'Within 1 Standard Deviation')
    plt.plot(js,np.round(np.array(perc_wihtin_acc)*100,2), linewidth = 1, c= 'b', label = 'Within 0.05 lDDT')
    plt.plot(js,[68.27]*len(js),'--', linewidth = 1, c= 'lightseagreen', label = 'Normal Distribution')
    plt.xlabel(xlabel)
    plt.ylabel('%')
    plt.legend()
    plt.ylim([0,100])
    plt.xticks([0,10,20,30,40,50,60,70,80,90,100])
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    pdb.set_trace()

    #Plot gradients
    fig = plt.figure(figsize=(9/2.54,9/2.54)) #set figsize
    plt.scatter(js, gradients,s=20)
    plt.plot(js, gradients, linewidth = 3)
    plt.ylabel('gradient')
    #plt.ylim(grad_ylims[score])
    plt.xlim([0,9.1])
    plt.xticks([0,1,2,3,4,5,6,7,8,9])
    if score == 'lddt_scores':
        plt.ylim([-0.025, 0.025])
    plt.xlabel(xlabel)
    fig.savefig(outdir+'gradient_running_'+suffix, format = 'svg')
    #Plot Point distribution
    fig = plt.figure(figsize=(9/2.54,9/2.54)) #set figsize
    plt.plot(js, perc_points, linewidth = 4)
    plt.xlabel(xlabel)
    plt.ylabel('% of points')
    #plt.xlim([2,10])
    #plt.xticks([2,3,4,5,6,7,8,9,10])
    fig.savefig(outdir+'perc_points_running_'+suffix, format = 'svg')

    av_df['ML '+cardinality[1:]+' distance'] = js
    av_df[score+aln_type] = avs
    av_df[score+aln_type+'_std'] = stds
    return av_df


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
get_one = bool(args.get_one[0])

aln_types = ['_straln','_seqaln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0], 'DIFFC':[0,1], 'TMscore': [0.2,1.0]}

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


#rename TMscore cols
hgroupdf = hgroupdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})
topdf = topdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})

cardinality = '_AA20'
av_df = pd.DataFrame()
for score in ['lddt_scores', 'TMscore', 'DIFFC', 'RMSD','DIFFSS', 'DIFF_ACC']:
    for aln_type in aln_types:
        ylim = ylims[score]
        av_df = ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir, av_df)
av_df.to_csv(outdir+'av_df.csv')
