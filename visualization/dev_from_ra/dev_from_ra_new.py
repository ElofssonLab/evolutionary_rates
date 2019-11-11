#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy import stats

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

parser.add_argument('--avdf', nargs=1, type= str,
default=sys.stdin, help = 'Dataframe with averages.')

def dev_from_av(avdf, df, score, aln_type):
    '''Calculate avearge dev from total running average within group and significance
    '''

    if cardinality == '_AA20':
            cardinality = ''

    avs = [] #Save average score
    avdevs = [] #Save deviations
    step = 0.1

    total_dists = avdf['ML '+cardinality[1:]+' distance']
    total_scores = avdf[score+aln_type]

    mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
    start = np.round(min(mldists),2)
    pdb.set_trace()
    end = min(max(mldists), 6) #End at 6
    scores = np.asarray(df[score+aln_type])

    for j in np.arange(min(mldists)+step,end+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)

        tav = total_avs[np.round(j-step/2, 1)] #total average in current interval
        avdevs.append(av-tav)


    #Do a t-test to calculate if the deviation is significant for each H-group
    #The variances will be unequal, as each H-group only encompasses very feq points
    #True values = 0, no deviation from total average
    truevals = np.zeros(len(df)) #size should not matter since zeros
    statistic, pvalue = stats.ttest_ind(avdevs, truevals, equal_var = False)


    return np.average(avdevs), pvalue


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
avdf = pd.read_csv(args.avdf[0])

cardinality = '_AA20'

#rename col
hgroupdf = hgroupdf.rename(columns={'H_group':'group'})
df = pd.concat([topdf, hgroupdf])
pdb.set_trace()
for score in ['RMSD','DIFFSS', 'DIFF_ACC', 'lddt_scores']:
    for aln_type in aln_types:

        av_from_line, pvalue = dev_from_av(avdf, df, score, aln_type)
