#! /usr/bin/env python3
# -*- coding: utf-8 -*-


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

parser.add_argument('--df', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with pdb files.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')

parser.add_argument('--plot_percentage', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot percentages or not.')

def runnning_average(outdir, complete_df, aln_type, score, cardinality, plot_gradients, plot_percentage, plot_num, pdf, fig,):
    '''Produce running average plots for df
    '''

    plt.rc('axes', titlesize=10, labelsize=10) #set title and label text sizes
    xlabel = 'MLAAdist'+cardinality+aln_type
    classes = {1:'Alpha', 2: 'Beta', 3: 'Alpha Beta', 4: 'Few 2ndary structures', 'total': 'Total'}
    colors = {1: 'royalblue', 2: 'k', 3: 'green', 4: 'violet', 'total': 'r'}
    sizes = {}
    sizes['total'] = 100

    if cardinality == '_AA20':
        cardinality = ''
    #Plot total average
    js = {} #Save dists
    perc_points = {}
    total_avs = {}
    step = 0.1 #what they used 2009
    df = complete_df
    perc_points['total'] = []
    gradients = {} #Save gradients

    for uca in [1.,2.,3.,4.]:
        perc_points[uca] = []
        js[uca] = []
        df = complete_df[complete_df['C._x']==uca]
        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score+aln_type])
        #plt.scatter(mldists, scores, color = 'wheat') #add all points
        avs = [] #Save average score

        for j in np.arange(min(mldists)+step,max(mldists)+step,step):
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
            below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step] #all below j, but at least j-step
            cut_scores = np.asarray(below_df[score+aln_type])
            if cut_scores.size == 0: #If no scores
                continue
            av= np.average(cut_scores)
            avs.append(av)
            js[uca].append(j-step/2) #Should be in the middle of interval

            #Save av distances and perc points
            perc_points[uca].append(len(below_df)/len(complete_df)*100)

        perc = np.round(len(df)*100/len(complete_df),2)
        sizes[uca] = perc

        #Include derivatives
        grads = np.gradient(avs)
        gradients[uca] = grads

        plt.plot(js[uca], avs, color =colors[int(uca)], linewidth = 1)


    #Meta for ra
    plt.ylabel('Running average '+ score)
    plt.title('Running average plot '+cardinality[1:]+' '+aln_type[1:])
    plt.xlim([0,6])
    plt.ylim([0,4])


    #Plot gradients
    if plot_gradients == True:
        plot_num+=3
        plt.subplot(plot_num) #set plot_num
        for i in [1.,2.,3.,4.]:
            plt.scatter(js[i], gradients[i],s=1, color =colors[int(i)] )
            plt.ylabel('Gradients')

        plt.xlim([0,6])
        plt.ylim([-0.1,0.1])


    if plot_percentage == True:
        plot_num+=3
        plt.subplot(plot_num) #set plot_num
        for label in [1.,2.,3.,4.]:
            plt.plot(js[label], perc_points[label], label = classes[label]+', '+str(sizes[label])+' %', color = colors[label], linewidth = 1)
        plt.xlabel(xlabel)
        plt.ylabel('% of points')
        plt.legend(loc = 'best')
        plt.xlim([0,6])
        plt.ylim([0,2])


    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
plot_gradients = args.plot_gradients[0]
plot_percentage = args.plot_percentage[0]
score = 'RMSD'
cardinality = '_AA20'
#Define pdf
pdf = PdfPages(outdir+score+cardinality+'.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize
plot_num = 331 #3*3 = 9 plots

for aln_type in ['_seqaln', '_straln']:
    pdf, fig = runnning_average(outdir, df, aln_type, score, cardinality, plot_gradients, plot_percentage, plot_num, pdf, fig)
    plot_num += 2

pdf.savefig(fig)
pdf.close()
