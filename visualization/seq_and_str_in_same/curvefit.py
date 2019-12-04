#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import Counter

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots a running average and its curvefit.''')

parser.add_argument('--avdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--avdf1', nargs=1, type= str,
default=sys.stdin, help = 'path to avdf with one pair per H-group from dataset 1.')

parser.add_argument('--topdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')


###FUNCTIONS###
def plot_poly(df):
    plots
#####MAIN#####
args = parser.parse_args()
avdf = pd.read_csv(args.avdf[0])
avdf1 = pd.read_csv(args.avdf1[0])
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
#concat dfs
catdf = pd.concat([topdf, hgroupdf])

x = np.array(avdf['ML  distance'])
y = np.array(avdf['lddt_scores_straln'])
x1 = np.array(avdf1['ML  distance'])
y1 = np.array(avdf1['lddt_scores_straln'])
#Fit polyline
z = np.polyfit(x, y, deg = 3)
p = np.poly1d(z)
z1 = np.polyfit(x1, y1, deg = 3)
p1 = np.poly1d(z1)



#Get onepairs
#set random seed
np.random.seed(42)
#get one pair per H-group from hgroupdf
groups = [*Counter(hgroupdf['group']).keys()]
one_pair_df = pd.DataFrame(columns = hgroupdf.columns)
for g in groups:
    partial_df = hgroupdf[hgroupdf['group']==g]
    i = np.random.randint(len(partial_df), size = 1)
    start =  partial_df.index[0]
    selection = partial_df.loc[start+i]
    one_pair_df = one_pair_df.append(selection)
#concat dfs
catdf1 = pd.concat([topdf,one_pair_df])

#Plot
matplotlib.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(10,10)) #set figsize
plt.scatter(catdf['MLAAdist_straln'], catdf['lddt_scores_straln'], label = 'Dataset 4', s= 1, c = 'b', alpha = 0.2)
plt.scatter(catdf1['MLAAdist_straln'], catdf1['lddt_scores_straln'], label = 'Dataset 5', s= 1, c = 'r', alpha = 0.2)
plt.plot(x,y, label = 'Running average Dataset 4',linewidth = 3, c= 'b')
plt.plot(x,p(x), label = '3 dg polynomial fit Dataset 4',linewidth = 3, c= 'deepskyblue')
plt.plot(x1,y1, label = 'Running average Dataset 5',linewidth = 3, c= 'r')
plt.plot(x1,p1(x1), label = '3 dg polynomial fit Dataset 5',linewidth = 3, c= 'mediumvioletred')

plt.legend(markerscale=10)
plt.ylim([0.2,1])
plt.xlim([0,9.1])
plt.xticks([0,1,2,3,4,5,6,7,8,9])
plt.xlabel('ML AA20 distance')
plt.ylabel('lDDT score')

fig.savefig(outdir+'curvefit.png', format = 'png')
print('Dataset 4',p)
print('Dataset 5',p1)

#Assess error towards polynomial
e=np.average(np.absolute(p(np.array(catdf['MLAAdist_straln']))-np.array(catdf['lddt_scores_straln'])))
pdb.set_trace()
print('Average error Dataset 4:', e)
e1=np.average(np.absolute(p1(np.array(catdf1['MLAAdist_straln']))-np.array(catdf1['lddt_scores_straln'])))
print('Average error Dataset 5:', e1)
