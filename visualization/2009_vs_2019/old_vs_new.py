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

parser.add_argument('--realigndf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')

parser.add_argument('--get_one', nargs=1, type= int,
default=sys.stdin, help = 'Get one pair from each H-group (1) or all (0).')


def ra_different(catdf, epsdf, realigndf, aln_type, score, cardinality, calc, ylim, outdir):
    '''Produce running average plots for df using 0.1 seqdist interval
    '''

    #Double column = 183 mm, maxfont = 8
    matplotlib.rcParams.update({'font.size': 7})
    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    xlabel = 'ML '+cardinality[1:]+' distance'
    grad_ylims = {'RMSD':[-0.15,0.15]}


    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''

    cat_avs = [] #Save average score
    cat_stds = [] #Save std deviation
    eps_avs = [] #Save average score
    eps_stds = [] #Save std deviation
    realign_avs = []


    js = [] #Save dists


    cat_mldists = np.asarray(catdf['MLAAdist'+cardinality+aln_type])
    cat_scores = np.asarray(catdf[score+aln_type])
    eps_mldists = np.asarray(epsdf['MLAAdist'+cardinality+aln_type])
    eps_scores = np.asarray(epsdf[score+aln_type])

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

        below_realigndf =realigndf[realigndf['MLAAdist']<j]
        below_realigndf = below_realigndf[below_realigndf['MLAAdist']>=j-step]
        cut_realign_scores = np.asarray(below_realigndf[score])

        if calc == 'average':
            catav= np.average(cut_cat_scores)
            epsav= np.average(cut_eps_scores)
            realignav= np.average(cut_realign_scores)
        if calc == 'median':
            catav= np.median(cut_cat_scores)
            epsav= np.median(cut_eps_scores)
            realignav= np.median(cut_realign_scores)

        cat_avs.append(catav)
        eps_avs.append(epsav)
        realign_avs.append(realignav)

        cat_std = np.std(cut_cat_scores)
        eps_std = np.std(cut_eps_scores)

        cat_stds.append(cat_std)
        eps_stds.append(eps_std)

        js.append(np.round(j-step/2,2))

    #Scatter
    ax.scatter(cat_mldists, cat_scores, c = 'cornflowerblue', s=0.1, linewidth = 2, label = 'Broad Dataset', alpha = 0.2)
    ax.scatter(eps_mldists, eps_scores, c = 'lightcoral',  s=0.1, linewidth = 2, label = '2009', alpha = 0.5)
    #Plot RA
    ax.plot(js, cat_avs, c = 'b', linewidth = 2, label = 'Broad Dataset')
    ax.plot(js, eps_avs, c = 'r', linewidth = 2, label = '2009')
    ax.plot(js, realign_avs, c = 'purple', linewidth = 2, label = '2009 TMalign')

    ax.set_ylabel(score)
    ax.legend(markerscale=5,fancybox=True, framealpha=0.5)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_xlim([0,6.1])
    ax.set_xticks([0,1,2,3,4,5,6])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'running_'+suffix, format = 'svg')
    fig.savefig(outdir+'running_'+suffix+'.png', format = 'png')
    plt.close()


    #Plot gradients
    fig, ax = plt.subplots(figsize=(9/2.54,9/2.54))
    gradients = np.gradient(eps_avs)
    ax.plot(js, gradients, linewidth = 1, c='r', label = '2009')
    smoothed_grads = ndimage.gaussian_filter1d(gradients, 2)
    ax.plot(js, smoothed_grads, '--', label = 'Gaussian KDE',linewidth = 2, c= 'r') #Plot gradients of gaussian kde

    gradients = np.gradient(cat_avs)
    ax.plot(js, gradients, linewidth = 1, c='b', label = 'Broad Dataset')
    smoothed_grads = ndimage.gaussian_filter1d(gradients, 2)
    ax.plot(js, smoothed_grads, '--', label = 'Gaussian KDE',linewidth = 2, c= 'b') #Plot gradients of gaussian kde

    gradients = np.gradient(realign_avs)
    ax.plot(js, gradients, linewidth = 1, c='purple', label = '2009 TMalign')
    smoothed_grads = ndimage.gaussian_filter1d(gradients, 2)
    ax.plot(js, smoothed_grads, '--', label = 'Gaussian KDE',linewidth = 2, c= 'purple') #Plot gradients of gaussian kde

    ax.set_ylabel('gradient')
    #plt.ylim(grad_ylims[score])
    ax.set_xlim([0,6.1])
    ax.set_xticks([0,1,2,3,4,5,6])
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
epsdf = pd.read_csv(args.epsdf[0])
realigndf = pd.read_csv(args.realigndf[0])
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
        ra_different(catdf, epsdf, realigndf, aln_type, score, cardinality, calc, ylim, outdir+score+aln_type+'/')
