#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import os
import argparse
from scipy import ndimage

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

def ridge_plot(df):
    '''Create a ridgeplot in pd
    'ML distance', 'Average', 'Interval scores'
    '''
    #fig = plt.figure(figsize=(9/2.54,9/2.54))
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'figure.figsize':(9/2.54,9/2.54)})
    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row="ML distance", hue="ML distance", aspect=15, height=0.5, palette=pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "lDDT score", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    g.map(sns.kdeplot, "lDDT score", clip_on=False, color="w", lw=2, bw=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)


    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)


    g.map(label, "lDDT score")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    pdb.set_trace()

    return None

def ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir, av_df):
    '''Produce running average plots for df using 100 steps with an equal amount of points in each
    '''

    #Double column = 183 mm, maxfont = 8
    matplotlib.rcParams.update({'font.size': 7})
    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025], 'DIFFC':[-0.04, 0.04], 'TMscore':[-0.025, 0.025]}


    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''
    avs = [] #Save average score
    stds = [] #Save std deviation
    js = [] #Save dists
    scores_per_interval = [] #Save scores per interval
    js_per_interval = []
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
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54)) #set figsize


    #Sort df by x-value
    df = df.sort_values(by=['MLAAdist'+cardinality+aln_type], ascending=True)
    step = 0.1

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
        perc_points.append(len(below_df)/len(df))
        #Get % points within 0.05
        outside = np.where(cut_scores<av-0.05)[0].size+np.where(cut_scores>av+0.05)[0].size
        if cut_scores.size ==0: #If no points in interval
            continue
        perc_within_acc.append((cut_scores.size-outside)/cut_scores.size)
        #Get % points within std
        outside_std = np.where(cut_scores<av-std)[0].size+np.where(cut_scores>av+std)[0].size
        perc_within_std.append((cut_scores.size-outside_std)/cut_scores.size)
        #Save scores in interval
        if j%1 == 0 or j == 0.1:
            scores_per_interval.extend(cut_scores)
            js_per_interval.extend([js[-1]]*cut_scores.size)
    #Include derivatives
    gradients = np.gradient(avs)
    #Fit polyline
    z = np.polyfit(js, avs, deg = 3)
    p = np.poly1d(z)
    #Plot RA
    l1 = ax.plot(js, avs, linewidth = 2, c = 'g', label = 'Running average')
    #ax.scatter(hgroup_mldists, hgroup_scores, s = 0.1, c = 'lightseagreen', alpha = 0.5, label = '95% Dataset')
    #ax.scatter(top_mldists, top_scores, s = 0.1, c = 'royalblue', alpha = 1.0, label = 'Topology Dataset')
    sns.kdeplot(top_mldists, top_scores,  shade=True, shade_lowest = False, cmap = 'Blues')
    sns.kdeplot(hgroup_mldists, hgroup_scores,  shade=True, shade_lowest = False, cmap = 'Greens')
    #plot stddev
    ax.plot(js, np.array(avs)+np.array(stds), '--', c = 'g', linewidth = 1) #positive stds
    l2 = ax.plot(js, np.array(avs)-np.array(stds), '--', c = 'g', linewidth = 1, label = 'Standard deviation') #negative stds
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score')
    else:
        ax.set_ylabel(score)
    #Custom legend
    patch1 = mpatches.Patch(color='lightseagreen', label='95% dataset')
    patch2 = mpatches.Patch(color='royalblue', label='Topology dataset')
    ax.legend(handles=[patch1, patch2, l1[0], l2[0]],markerscale=5, frameon=False)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'running_'+suffix, format = 'svg')
    fig.savefig(outdir+'running_'+suffix+'.png', format = 'png')
    plt.close()

    #Plot spread within each interval in ridgelpot
    #spread_df = pd.DataFrame(list(zip(js_per_interval, scores_per_interval)), columns = ['ML distance', 'lDDT score'])
    #ridge_plot(spread_df)

    #Plot polynomial fit
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    #ax.scatter(hgroup_mldists, hgroup_scores, s = 0.1, c = 'lightseagreen', alpha = 0.8, label = '95% Dataset')
    #ax.scatter(top_mldists, top_scores, s = 0.1, c = 'royalblue', alpha = 1.0, label = 'Topology Dataset')
    sns.kdeplot(top_mldists, top_scores,  shade=True, shade_lowest = False, cmap = 'Blues')
    sns.kdeplot(hgroup_mldists, hgroup_scores,  shade=True, shade_lowest = False, cmap = 'Greens')

    ax.plot(js, avs, linewidth = 2, c = 'g', label = 'Running average')
    ax.plot(js,p(js), label = '3 dg polynomial fit',linewidth = 1, c= 'b')
    plt.title('Broad Dataset')
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score')
    else:
        ax.set_ylabel(score)
    ax.legend(markerscale=7, fancybox=True, frameon=False)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'polynomial_'+suffix, format = 'svg')
    fig.savefig(outdir+'polynomial_'+suffix+'.png', format = 'png')
    plt.close()

    #Plot histograms of x and y vals
    #str score
    fig, ax = plt.subplots(figsize=(3/2.54,9/2.54))
    ax.hist(topdf[score+aln_type], color = 'royalblue', orientation = 'horizontal', density = True, alpha = 0.8, bins = 30)
    ax.hist(hgroupdf[score+aln_type], color = 'lightseagreen', orientation = 'horizontal', density = True, alpha = 0.5, bins = 30)
    ax.set_xlabel(xlabel)
    ax.set_ylim(ylim)
    ax.invert_xaxis()
    ax.axis('off')
    #fig.tight_layout()
    fig.savefig(outdir+'strhist_'+suffix, format = 'svg')
    fig.savefig(outdir+'strhist_'+suffix+'.png', format = 'png')
    plt.close()
    #seqdist
    fig, ax = plt.subplots(figsize=(9/2.54,3/2.54))
    ax.hist(topdf['MLAAdist'+cardinality+aln_type], color = 'royalblue', density = True, alpha = 0.8, bins = 30)
    ax.hist(hgroupdf['MLAAdist'+cardinality+aln_type], color = 'lightseagreen', density = True, alpha = 0.5, bins = 30)
    ax.set_xlabel(xlabel)
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.invert_yaxis()
    ax.axis('off')
    #fig.tight_layout()
    fig.savefig(outdir+'seqhist_'+suffix, format = 'svg')
    fig.savefig(outdir+'seqhist_'+suffix+'.png', format = 'png')
    plt.close()
    #Plot Percent within accuracy and std
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    ax.plot(js,np.round(np.array(perc_within_std)*100,2), linewidth = 1, c= 'g',  label = 'Within 1 Std')
    ax.plot(js,np.round(np.array(perc_within_acc)*100,2), linewidth = 1, c= 'b', label = 'Within 0.05 lDDT')
    ax.plot(js,[68.27]*len(js),'--', linewidth = 1, c= 'lightseagreen')
    ax.set_ylabel('% of Pairs')
    ax.legend(frameon=False)
    ax.set_ylim([20,100])
    ax.set_yticks([20,30,40,50,60,70,80,90,100])
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_xlabel(xlabel)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'within_dev_'+suffix, format = 'svg')
    fig.savefig(outdir+'within_dev_'+suffix+'.png', format = 'png')
    plt.close()


    #Plot std against sequence distnace
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    ax.plot(js,stds, linewidth = 1, c= 'g', label = 'Standard Deviation')
    ax.plot(js,[0.05]*len(js),'--', linewidth = 1, c= 'b')
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score')
    else:
        ax.set_ylabel(score)
    ax.legend(frameon=False)
    ax.set_ylim([0,0.16])
    ax.set_yticks(np.arange(0,0.17,0.02))
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_xlabel(xlabel)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'seqdist_std'+suffix, format = 'svg')
    fig.savefig(outdir+'seqdist_std'+suffix+'.png', format = 'png')
    plt.close()

    #Plot relative std (std/mean) against sequence distnace
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    ax.plot(js,np.array(stds)/np.array(avs), linewidth = 1, c= 'g', label = 'Relative Standard Deviation')
    ax.plot(js,[0.05]*len(js),'--', linewidth = 1, c= 'b')
    ax.set_ylabel('Std/Mean')
    ax.legend(frameon=False)
    ax.set_ylim([0,0.24])
    ax.set_yticks(np.arange(0,0.25,0.02))
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_xlabel(xlabel)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'seqdist_rel_std'+suffix, format = 'svg')
    fig.savefig(outdir+'seqdist_rel_std'+suffix+'.png', format = 'png')
    plt.close()

    #Plot gradients
    fig, ax = plt.subplots(figsize=(9/2.54,6/2.54))
    #ax.scatter(js, gradients,s=2)
    ax.plot(js, gradients, linewidth = 1)
    smoothed_grads = ndimage.gaussian_filter1d(gradients, 2)
    ax.plot(js, smoothed_grads, label = '1D Gaussian KDE',linewidth = 1, c= 'indigo') #Plot gradients of polyline
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score gradient')
    else:
        ax.set_ylabel(score+' gradient')
    #plt.ylim(grad_ylims[score])
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_ylim(grad_ylims[score])
    ax.set_xlabel(xlabel)
    ax.legend(frameon=False)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'gradient_running_'+suffix, format = 'svg')
    fig.savefig(outdir+'gradient_running_'+suffix+'.png', format = 'png')
    plt.close()

    #Plot Point distribution
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54)) #set figsize
    ax.plot(js, perc_points, linewidth = 2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('% of points')
    #plt.xlim([2,10])
    #plt.xticks([2,3,4,5,6,7,8,9,10])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig(outdir+'perc_points_running_'+suffix, format = 'svg')
    fig.savefig(outdir+'perc_points_running_'+suffix+'.png', format = 'png')

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
ylims = {'RMSD':[0,5], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0], 'DIFFC':[0,1], 'TMscore': [0.2,1.0]}

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
        try:
            os.mkdir(outdir+score+aln_type)
        except:
            print('Directory exists')
        ylim = ylims[score]
        av_df = ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir+score+aln_type+'/', av_df)
av_df.to_csv(outdir+'av_df.csv')
