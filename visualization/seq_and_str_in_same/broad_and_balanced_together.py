#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sys
import os
import argparse

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots running averages and 3dg polynomial fits.''')

parser.add_argument('--avdf_all', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--avdf_one', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

#####MAIN#####
args = parser.parse_args()
avdf_all= pd.read_csv(args.avdf_all[0])
avdf_one= pd.read_csv(args.avdf_one[0])
outdir = args.outdir[0]


#Plot
matplotlib.rcParams.update({'font.size': 7})
fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
ax.plot(avdf_all['ML  distance'], avdf_all['lddt_scores_straln'], linewidth = 2, c = 'cornflowerblue', label = 'Broad')
ax.plot(avdf_one['ML  distance'], avdf_one['lddt_scores_straln'], linewidth = 2, c = 'mediumpurple', label = 'Balanced Broad')

#Polynomials
z = np.polyfit(avdf_all['ML  distance'], avdf_all['lddt_scores_straln'], deg = 3)
p = np.poly1d(z)
ax.plot(avdf_all['ML  distance'],p(avdf_all['ML  distance']),linewidth = 1, c= 'mediumblue')
z = np.polyfit(avdf_one['ML  distance'], avdf_one['lddt_scores_straln'], deg = 3)
p = np.poly1d(z)
ax.plot(avdf_one['ML  distance'],p(avdf_one['ML  distance']),linewidth = 1, c= 'rebeccapurple')

#Layout
ax.set_xlim([0,9.1])
ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
ax.set_ylim([0.2,1])
ax.set_yticks(np.arange(0.2,1.1,0.1))
xlabel='ML AA20 distance'
ylabel='lDDT score'
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.title('Broad and Balanced Broad')
ax.legend(markerscale=5,fancybox=True, framealpha=0.5)
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout()
#save
fig.savefig(outdir+'broad_and_balanced.svg', format = 'svg')
fig.savefig(outdir+'broad_and_balanced.png', format = 'png')
plt.close()
