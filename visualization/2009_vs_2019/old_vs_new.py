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

parser.add_argument('--epsdf', nargs=1, type= str,
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

def ra_different(catdf, epsdf, aln_type, score, cardinality, calc, ylim, outdir+score+aln_type+'/'):
    '''Produce running average plots for df using 0.1 seqdist interval
    '''

    #Double column = 183 mm, maxfont = 8
    matplotlib.rcParams.update({'font.size': 7})
    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.1,0.1]}


    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''

    cat_avs = [] #Save average score
    cat_stds = [] #Save std deviation
    cat_avs = [] #Save average score
    cat_stds = [] #Save std deviation


    js = [] #Save dists


    cat_mldists = np.asarray(catdf['MLAAdist'+cardinality+aln_type])
    cat_scores = np.asarray(catdf[score+aln_type])
    eps_mldists = np.asarray(epsdf['MLAAdist'+cardinality+aln_type])
    hgroup_scores = np.asarray(epsdf[score+aln_type])

    mldists = np.append(cat_mldists, eps_mldists)
    scores = np.append(cat_scores, eps_scores)
    step = 0.1

    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54)) #set figsize



    for j in np.arange(min(mldists)+step,6+step,step):
        below_epsdf = epsdf[epsdf['MLAAdist'+cardinality+aln_type]<j]
        below_epsdf = below_epsdf[below_epsdf['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_eps_scores = np.asarray(below_epsdf[score+aln_type])

        below_catdf =catdf[catdf['MLAAdist'+cardinality+aln_type]<j]
        below_catdf = below_catdf[below_catdf['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_cat_scores = np.asarray(below_catdf[score+aln_type])
        if calc == 'average':
            catav= np.average(cut_cat_scores)
            epsav= np.average(cut_eps_scores)
        if calc == 'median':
            catav= np.median(cut_cat_scores)
            epsav= np.median(cut_eps_scores)

        cat_avs.append(catav)
        eps_avs.append(apsav)

        cat_std = np.std(cut_cat_scores)
        eps_std = np.std(cut_eps_scores)
        cat_stds.append(cat_std)
        eps_stds.append(eps_std)

        js.append(np.round(j-step/2,2))

    #Include derivatives
    gradients = np.gradient(eps_avs)

    #Plot RA
    ax.plot(js, avs, linewidth = 2, c = 'g', label = 'Running average')
    #plot stddev
    ax.plot(js, eps_avs, '--', c = 'g', linewidth = 1)
    ax.plot(js, cat_avs, '--', c = 'g', linewidth = 1)
    ax.set_ylabel(score)
    ax.legend(markerscale=5,fancybox=True, framealpha=0.5)
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


    #Plot gradients
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    #ax.scatter(js, gradients,s=2)
    ax.plot(js, gradients, linewidth = 1)
    smoothed_grads = ndimage.gaussian_filter1d(gradients, 2)
    ax.plot(js, smoothed_grads, label = '1D Gaussian KDE',linewidth = 1, c= 'indigo') #Plot gradients of polyline
    ax.set_ylabel('gradient')
    #plt.ylim(grad_ylims[score])
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_ylim(grad_ylims[score])
    ax.set_xlabel(xlabel)
    ax.legend()
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'gradient_running_'+suffix, format = 'svg')
    fig.savefig(outdir+'gradient_running_'+suffix+'.png', format = 'png')
    plt.close()


    return None


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
epsdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
get_one = bool(args.get_one[0])

aln_types = ['_straln']
ylims = {'RMSD':[0,4]}

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

catdf = pd.concat([topdf, hgroupdf])
#rename
epsdf = epsdf.rename(columns={'seqdist':'MLAAdist_straln', 'RMSD':'RMSD_straln'})
cardinality = '_AA20'
av_df = pd.DataFrame()
for score in ['RMSD']:
    for aln_type in aln_types:
        try:
            os.mkdir(outdir+score+aln_type)
        except:
            print('Directory exists')
        ylim = ylims[score]
        ra_different(catdf, epsdf, aln_type, score, cardinality, calc, ylim, outdir+score+aln_type+'/')
