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

parser.add_argument('--hgroupdf95', nargs=1, type= str,
default=sys.stdin, help = 'path to 95df.')

parser.add_argument('--hgroupdf20', nargs=1, type= str,
default=sys.stdin, help = 'path to 20df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')











#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf95 = pd.read_csv(args.hgroupdf95[0])
hgroupdf20 = pd.read_csv(args.hgroupdf20[0])
outdir = args.outdir[0]

dfs = [topdf, hgroupdf95, hgroupdf20]
suffix = ['_topdf.svg', '_hgroupdf95.svg', '_hgroupdf20.svg']

labels = {1.:'Mainly Alpha', 2.: 'Mainly Beta', 3.:'Alpha Beta', 4.: 'Few SS'}
colors = {1.: 'royalblue', 2.: 'k', 3.: 'green', 4.: 'violet'}
matplotlib.rcParams.update({'font.size': 22})
for i in range(3):
    df = dfs[i]
    print(len(df))

    for aln_type in ['_seqaln', '_straln']:
        fig = plt.figure(figsize=(10,10)) #set figsize
        for key in labels:
            class_df = df[df['C._x']==key]
            label = labels[key]
            sns.distplot(class_df['aln_len'+aln_type], color = colors[key], hist = False, label = label, linewidth = 2)
            plt.xlabel('Aligned length')
            plt.ylabel('Density')
            plt.xlim([0,300])
            plt.ylim([0,0.04])
            plt.legend()
        fig.savefig(outdir+'kde'+aln_type+suffix[i], format = 'svg')
        plt.close()


pdb.set_trace()
