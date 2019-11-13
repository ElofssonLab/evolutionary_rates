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
suffix = ['topdf.svg', 'hgroupdf95.svg', 'hgroupdf20.svg']
for i in range(3):
    df = dfs[i]
    print(len(df))
    fig = plt.figure(figsize=(10,10)) #set figsize
    ax = sns.pairplot(data = df, vars = ['aln_len_seqaln', 'aln_len_straln'], hue = 'C._x')
    ax.savefig(outdir+'pairplot'+suffix[i], format = 'svg')


pdb.set_trace()
