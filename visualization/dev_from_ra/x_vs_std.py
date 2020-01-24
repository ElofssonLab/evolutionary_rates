#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy import stats
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols
import time
'''
When using statsmodels in scientific publication, please consider using the following citation:
Seabold, Skipper, and Josef Perktold. “Statsmodels: Econometric and statistical modeling with python.”
Proceedings of the 9th Python in Science Conference. 2010.
'''

#Custom
from percent_ss import parse_ss
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that calculates deviations for running averages.''')

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

def AA6_distribution(df, aln_type):
    '''Calculate AA6 distribution in sequence (gaps have their own state as well)
    Groupings: K=[KR],D=[DE],Y=[YWFH],T=[TSQN],C=[CVMLIA], P=[PG], "-"="-"
    '''
    #Save percentages:
    percentage_dict = {
    'K1'+aln_type:[],
    'K2'+aln_type:[],
    'D1'+aln_type:[],
    'D2'+aln_type:[],
    'Y1'+aln_type:[],
    'Y2'+aln_type:[],
    'T1'+aln_type:[],
    'T2'+aln_type:[],
    'C1'+aln_type:[],
    'C2'+aln_type:[],
    'P1'+aln_type:[],
    'P2'+aln_type:[],
    '-1'+aln_type:[],
    '-2'+aln_type:[]
    }


    #reset index
    df = df.reset_index()
    t1 = time.time()
    for i in range(len(df)):
        row = df.iloc[i]
        l1 = len(row['AA6_1_straln'])
        l2 = len(row['AA6_2_straln'])
        s1 = Counter(row['AA6_1_straln'])
        s2 = Counter(row['AA6_2_straln'])

        #Get s1 ad s2 percentages
        for key in ['P', 'C', 'K', 'T', 'D', 'Y', '-']:
            try:
                percentage_dict[key+'1'+aln_type].append(s1[key]/l1)
            except:
                percentage_dict[key+'1'+aln_type].append(0)
            try:
                percentage_dict[key+'2'+aln_type].append(s2[key]/l2)
            except:
                percentage_dict[key+'2'+aln_type].append(0)

    print('Calculating AA6',time.time()-t1,'s')

    for key in percentage_dict:
        df[key] = percentage_dict[key]

    return df, percentage_dict.keys()

def dev_from_av(avdf, df, score, aln_type, cardinality):
    '''Calculate number of stds from total running average for each pair
    '''
    start = 0
    end = 6
    step = 0.1
    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_df = df[df['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]

        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        cut_scores = np.asarray(below_df[score+aln_type])

        #Get total average in interval
        x = np.round(j-step/2, 2)
        tav = avdf[avdf['ML  distance']==x][score+aln_type].values[0] #total average in current interval
        #Get std in interval
        std = avdf[avdf['ML  distance']==x][score+aln_type+'_std'].values[0]
        #Deviation
        dev = cut_scores-tav
        std_away = dev/std

    return np.average(avdevs), pvalue, statistic, js, avs


#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
avdf = pd.read_csv(args.avdf[0])

cardinality = '_AA20'

#Get topology from hgroupdf
tops = []
hgroups = [*hgroupdf['group']]
for hg in hgroups:
    hg = hg.split('.')
    tops.append(hg[0]+'.'+hg[1]+'.'+hg[2])

hgroupdf['C.A.T.'] = tops
#rename col
hgroupdf = hgroupdf.rename(columns={'group':'H_group'})
hgroupdf = hgroupdf.rename(columns={'C.A.T.':'group'})

#rename TMscore cols
hgroupdf = hgroupdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})
topdf = topdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})

catdf = pd.concat([topdf, hgroupdf])
#Rename class column
catdf = catdf.rename(columns={'C._x':'Class'})

#The ones should actually be zeros
catdf['RCO1']=catdf['RCO1'].replace([1], 0)
catdf['RCO2']=catdf['RCO2'].replace([1], 0)


#Save pvalues
top_metrics = pd.DataFrame()
top_metrics['Topology'] = topologies
cardinality = '' #AA20

for score in ['lddt_scores', 'TMscore', 'DIFFC', 'RMSD', 'DIFFSS', 'DIFF_ACC']:
    for aln_type in ['_straln', '_seqaln']:
        #select below 6 using seq or str
        catdf_s = catdf[catdf['MLAAdist'+aln_type]<=6]

        features = ['RCO', 'aln_len'+aln_type, 'l', 'percent_aligned'+aln_type,'P', 'C', 'K', 'T', 'D', 'Y',
        '-', 'L', 'S', 'H', 'CD'] #L,S,H = loop, sheet, helix

        catdf_s, perc_keys = AA6_distribution(catdf_s, aln_type) #Get AA6 frequencies
        catdf_s = parse_ss(catdf_s, aln_type) #Get % ss

        std = dev_from_av(avdf, df, score, aln_type, cardinality, 6)
