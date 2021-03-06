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
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
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



def ra_different(topdf, hgroupdf, aln_type, score, cardinality, calc, ylim, outdir):
    '''Produce running average plots for df
    '''

    suffix = calc+'_'+score+cardinality+aln_type+'_'+'.svg'
    colors = {'_seqaln': 'k', '_straln': 'r'}
    xlabel = cardinality[1:]+' ED'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025]}
    plt.rc('axes', titlesize=10) #set title and label text sizes

    gradients = {}

    #Plot total average for cardinality
    if cardinality == '_AA20':
        cardinality = ''
    avs = [] #Save average score
    js = [] #Save dists
    stds = [] #Save stds
    perc_points = []
    num_points = 0
    step = 0.1
    top_mldists = np.asarray(topdf['MLAAdist'+cardinality+aln_type])
    top_scores = np.asarray(topdf[score+aln_type])
    hgroup_mldists = np.asarray(hgroupdf['MLAAdist'+cardinality+aln_type])
    hgroup_scores = np.asarray(hgroupdf[score+aln_type])

    mldists = np.append(top_mldists, hgroup_mldists)
    scores = np.append(top_scores, hgroup_scores)
    df = pd.concat([topdf, hgroupdf])

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
        js.append(np.round(j-step/2,2))
        stds.append(np.std(cut_scores))
        perc_points.append(len(below_df)/len(df)*100)
        num_points += len(below_df[below_df['MLAAdist'+cardinality+aln_type]<=6])
    #Include derivatives
    gradients = np.gradient(avs)

    return js, avs, stds, perc_points, num_points, mldists, scores, gradients

def make_plots(results, cardinality, outdir, suffix):
    '''Make the plots'''
    matplotlib.rcParams.update({'font.size': 7})
    classes = {1.:'Mainly Alpha', 2.: 'Mainly Beta', 3.: 'Alpha Beta', 4.: 'Few SS'}
    colors = {1.: 'cornflowerblue', 2.: 'dimgray', 3.: 'mediumseagreen', 4.: 'darkorchid'}
    cmaps = {1.: 'Blues', 2.: 'Greys', 3.: 'Greens', 4.: 'Purples'}
    xlabel = cardinality[1:]+' ED'
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025]}

    for i in [1.,2.,3.,4.]:
        fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
        js, avs, stds, perc_points, num_points, mldists, scores, gradients = results[i]

        #ax.scatter3D(mldists, [i]*len(mldists),scores, c = colors[i], s = 0.05)

        sns.kdeplot(mldists, scores,  shade=True, shade_lowest = False, cmap = cmaps[i])
        #Plot RA
        ax.plot(js, avs, linewidth = 2, c = colors[i], label = 'Running average')
        ax.plot(js, np.array(avs)+np.array(stds), '--', linewidth = 1,  c = colors[i])
        ax.plot(js,np.array(avs)-np.array(stds), '--', linewidth = 1, c = colors[i], label = 'Standard deviation')
        #plt.scatter(mldists, scores, s = 0.05, c = colors[i], alpha = 0.7)
        plt.title(classes[i])
        ax.set_xlabel(xlabel)
        if score == 'lddt_scores':
            ax.set_ylabel('lDDT score')
        else:
            plt.ylabel(score)
        ax.set_ylim(ylim)
        ax.set_xlim([0,6.1])
        ax.set_xticks([0,1,2,3,4,5,6])
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(markerscale=7, frameon=False)
        fig.tight_layout()
        fig.savefig(outdir+'class_'+str(i)+'_running_'+suffix, format = 'svg')
        fig.savefig(outdir+'class_'+str(i)+'_running_'+suffix+'.png', format = 'png')
        plt.close()



    #Plot all together
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    for i in [1.,2.,3.,4.]:
        js, avs, stds, perc_points, num_points, mldists, scores, gradients = results[i]
        print(i, num_points)
        #Plot RA
        ax.plot(js, avs, linewidth = 2, c = colors[i], label = classes[i])

    ax.set_xlabel(xlabel)
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score')
    else:
        plt.ylabel(score)
    plt.title('All Classes')
    ax.set_ylim(ylim)
    ax.set_xlim([0,6.1])
    ax.set_xticks([0,1,2,3,4,5,6])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(markerscale=7, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir+'all_classes_running_'+suffix, format = 'svg')
    fig.savefig(outdir+'all_classes_running_'+suffix+'.png', format = 'png')
    plt.close()
    return None

def compare_classes(avdf, outfile):
    '''Compare the classes running averages
    '''
    classes = [1.,2.,3.,4.]
    f = open(outfile, 'w')
    f.write('Two-sided t-tests between ras of classes where both ras have values.\n')
    f.write('Alpha\tBeta\tAlphaBeta\tFewSS')
    for i in range(4):
        C1 = classes[i]
        f.write('\n')
        f.write('\t'*int(C1))
        for j in range(i+1,4):
            C2 = classes[j]
            test_df = pd.DataFrame() #Only want to compare points where both ras have values
            test_df[C1] = avdf[C1]
            test_df[C2] = avdf[C2]
            test_df = test_df.dropna() #Drop nans
            statistic, pvalue = stats.ttest_ind(test_df[C1], test_df[C2], equal_var = False)
            f.write(str(np.round(pvalue,3))+'\t')

    f.close() #Close file

#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
get_one = bool(args.get_one[0])

aln_types = ['_seqaln', '_straln']
ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0]}
cardinality = '_AA20'
results = {}

#set random seed
np.random.seed(42)
if get_one == True:
    #get one pair per H-group from hgroupdf
    groups = [*Counter(hgroupdf['H_group']).keys()]
    one_pair_df = pd.DataFrame(columns = hgroupdf.columns)
    for g in groups:
        partial_df = hgroupdf[hgroupdf['H_group']==g]
        i = np.random.randint(len(partial_df), size = 1)
        start =  partial_df.index[0]
        selection = partial_df.loc[start+i]
        one_pair_df = one_pair_df.append(selection)
    hgroupdf = one_pair_df

for score in ['lddt_scores']:#['RMSD','DIFFSS', 'DIFF_ACC', 'lddt_scores']:
    ylim = ylims[score]
    for aln_type in aln_types:
        avdf = pd.DataFrame() #Create dataframe to save avs
        for C in [1.,2.,3.,4.]: #Plot per class
            df1 = topdf[topdf['C._x']==C]
            df2 = hgroupdf[hgroupdf['C._x']==C]
            js, avs, stds, perc_points, num_points, mldists, scores, gradients = ra_different(df1, df2, aln_type, score, cardinality, calc, ylim, outdir)
            results[C] = [js, avs, stds, perc_points, num_points, mldists, scores, gradients]
            try:
             avdf[C] = avs
            except:
                pdb.set_trace()#Remember to have MLAAdist below 6 for stats calcs
        suffix = calc+'_'+score+cardinality+aln_type+'.svg'
        make_plots(results, cardinality, outdir, suffix)
        #Compare ras
        #outfile = outdir+score+cardinality+aln_type+'.pvals'
        #compare_classes(avdf, outfile)
