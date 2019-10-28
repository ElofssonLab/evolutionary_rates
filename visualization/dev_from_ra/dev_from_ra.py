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
from scipy import stats

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots the deviations from the running within H-groups.''')

parser.add_argument('--df', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either mean or average.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')

parser.add_argument('--H_groups', nargs=1, type= str,
                  default=sys.stdin, help = 'path to H_groups.')

parser.add_argument('--sequences', nargs=1, type= str,
                  default=sys.stdin, help = 'path to clustered sequences.')

def read_tsv(H_groups_file):
	'''Read ids and H-groups into dict
	'''

	H_groups = {} #Store H-groups and uids

	with open(H_groups_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1:]
			H_group = H_group[0]+'.'+H_group[1]+'.'+H_group[2]+'.'+H_group[3]
			H_groups[uid] = H_group


	return H_groups

def read_fasta(fasta_file):
	'''Read fasta file into dict
	'''

	sequences = {} #Store sequences
	with open(fasta_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			if line[0] == '>':
				uid = line.split('|')[2].split('/')[0]
			else:
				sequences[uid] = line

	return sequences

def original_length(sequences):
    '''Calculate the original sequence lengths
    '''
    lengths = {}
    for key in sequences:
        lengths[key] = len(sequences[key])

    return lengths
def get_groups(H_groups, sequences):
	'''Get H-group for each uid and group sequences accordingly
	'''


	grouped_sequences = {} #Sequences grouped by H-group

	for key in sequences:
		H_group = H_groups[key]
		sequence = sequences[key]

		if H_group not in grouped_sequences.keys(): #If not in new dict - add
			grouped_sequences[H_group] = [key + '/' + sequence]
		else:
			grouped_sequences[H_group].append(key + '/' + sequence) #Otherwise append

	return grouped_sequences

def count_groups():
    '''Get the original size for each H-group
    '''
    H_group_sizes = {}
    for uid in H_groups:
        group = H_groups[uid]
        if group not in H_group_sizes.keys():
            H_group_sizes[group] = 1
        else:
            H_group_sizes[group] += 1

def ra_different(df, aln_types, score, cardinality, calc):
    '''calculate running average for df
    '''

    ras= {} #Save running averages

    for i in range(2):
        aln_type = aln_types[i]
        #Plot total average for cardinality
        if cardinality == '_AA20':
            cardinality = ''
        avs = [] #Save average score
        js = [] #Save dists
        total_avs = {}
        step = 0.1
        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score+aln_type])

        for j in np.arange(min(mldists)+step,max(mldists)+step,step):
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
            below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
            cut_scores = np.asarray(below_df[score+aln_type])
            if calc == 'average':
                av= np.average(cut_scores)
            if calc == 'mean':
                av= np.mean(cut_scores)
            avs.append(av)
            js.append(j-step/2)
            total_avs[j-step] = av
        #Save running averages to dict
        ras[aln_type] = total_avs

    return ras

def distance_from_av(df, aln_type, score, cardinality, calc, tra, plot_num):
    '''Calculate the difference from the total running average for each H-group
    '''
    #Get unique groups
    H_groups = [*Counter(df['H_group']).keys()]
    #Set cardinality
    if cardinality == '_AA20':
        cardinality = ''
    #Set step
    step = 0.1

    #True values = 0, no deviation from total average
    truevals = np.zeros(105) #Maximum 105 pairs per H-group
    all_x = []
    all_y = []
    pvals = []
    mldists = [*tra.keys()] #Define thresholds as collected in ras
    for group in H_groups:
        hgroup_df = df[df['H_group']==group]
        scores = np.asarray(hgroup_df[score+aln_type])

        js = [] #save middle points in mldists
        avdevs = [] #save average deviations from line for each step

        #Only get the values where the H-group has points
        end = min(max(hgroup_df['MLAAdist'+cardinality+aln_type]), 6)

        for j in np.arange(min(mldists)+step, end+step, step):
            below_df = hgroup_df[hgroup_df['MLAAdist'+cardinality+aln_type]<j] #below j
            below_df = hgroup_df[hgroup_df['MLAAdist'+cardinality+aln_type]>=j-step] #above or equal to j-step
            cut_scores = np.asarray(hgroup_df[score+aln_type])


            tav = tra[j-step] #total average in current interval

            if calc == 'average':
                av= np.average(cut_scores)
            if calc == 'mean':
                av= np.mean(cut_scores)

            js.append(j-step/2)
            avdevs.append(av-tav)

        #Do a t-test to calculate if the deviation is significant for each H-group
        #The variances will be unequal, as each H-group only encompasses very feq points
        statistic, pvalue = stats.ttest_ind(avdevs, truevals, equal_var = False)
        pvals.append(pvalue)
        #Do a scatter plot of the deviation for each evdist
        plt.scatter(js, avdevs, s= 1)
        all_x.extend(js)
        all_y.extend(avdevs)

    plt.subplot(plot_num) #set plot_num
    plt.title(aln_type[1:])
    plt.ylabel('Average deviation from total average')
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])


    #New plot_num
    plot_num += 3
    plt.subplot(plot_num) #set plot_num
    sns.kdeplot(all_x, all_y, shade = True)
    plt.ylabel('Average deviation from total average')
    plt.xlabel('ML AA20')
    plt.xlim([0,6])
    plt.xticks([0,1,2,3,4,5,6])

    #New plot_num
    plot_num += 3
    plt.subplot(plot_num) #set plot_num
    plt.hist(pvals, bins = 30, log = True)
    plt.ylabel('Count')
    plt.xlabel('P-value')

    return all_x, all_y, pvals

#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
calc = args.calc[0]
plot_gradients = args.plot_gradients[0]
H_groups = read_tsv(args.H_groups[0]) #Read H-groups and uids
sequences = read_fasta(args.sequences[0]) #Read the clustered sequences
original_seqlens = original_length(sequences) #Get the original sequence lengths for the domains
grouped_sequences = get_groups(H_groups, sequences) #Group sequences


#Calculate total running averages
cardinality = '_AA20'
aln_types = ['_seqaln', '_straln']
score = 'RMSD'
ras = ra_different(df, aln_types, score, cardinality, calc)
#Calculate deviations from total ra and plot
plot_num = 331
for key in ras:
    fig = plt.figure(figsize=(10,10)) #set figsize
    tra = ras[key]
    aln_type = key
    all_x, all_y, pvals = distance_from_av(df, aln_type, score, cardinality, calc, tra, plot_num)

    #Save calculated values
    np.save(outdir+aln_type+'_all_x.npy', all_x)
    np.save(outdir+aln_type+'_all_y.npy', all_y)
    np.save(outdir+aln_type+'_pvals.npy', pvals)
    plot_num += 2

fig.savefig(outdir+score[1:]+'_'+calc+'_average_deviation.svg', format = 'svg')
pdb.set_trace()
