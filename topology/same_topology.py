#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import random
import os
import numpy as np

import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that filters entries in CATH
        according to a filter file containing PDB ids. All entries in each topology are then 
	randomly paired, using one entry per H-group and written to individual folders. 
	One folder per topology. Subfolders are created in batches of 200 entries within a topology.''')

parser.add_argument('--H_groups', nargs=1, type= str,
                  default=sys.stdin, help = 'path to H_groups.')
parser.add_argument('--sequences', nargs=1, type= str,
                  default=sys.stdin, help = 'path to clustered sequences.')
parser.add_argument('--failed_pdb_filter', nargs=1, type= str,
                  default=sys.stdin, help = 'path to uids that failed the pdb filter.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'path output directory.')
parser.add_argument('--logfile', nargs=1, type= str,
                  default=sys.stdin, help = 'logfile for stats.')


def read_tsv(H_groups_file):
	'''Read ids and H-groups into dict
	'''

	H_groups = {} #Store H-groups and uids
	topologies = {} #Store topologies and uids
	with open(H_groups_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1:]
			#Get topology
			topology = H_group[0]+'.'+H_group[1]+'.'+H_group[2]
			topologies[uid] = topology
			#Get H-group
			H_group = H_group[0]+'.'+H_group[1]+'.'+H_group[2]+'.'+H_group[3]
			H_groups[uid] = H_group

	return H_groups, topologies

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

def get_groups(topologies, sequences):
	'''Get topology for each uid and group sequences accordingly
	'''


	grouped_sequences = {} #Sequences grouped by topology

	for key in sequences:
		topology = topologies[key]
		sequence = sequences[key]

		if topology not in grouped_sequences.keys(): #If not in new dict - add
			grouped_sequences[topology] = [key + '/' + sequence]
		else:
			grouped_sequences[topology].append(key + '/' + sequence) #Otherwise append

	return grouped_sequences


def get_passed_uids(failed_pdb_filter, sequences, topologies, num_topologies_before):
    '''Get the uids that passed the pdb filter and group them in topologies
    '''

    #Get passed uids
    failed_uids = [*failed_pdb_filter[0]]
    passed_uids = []
    for uid in sequences:
        if uid in failed_uids:
            continue
        else:
            passed_uids.append(uid)

    #Group the passed uids by topology
    passed_uids_grouped = {} #Sequences grouped by topology
    for uid in passed_uids:
        topology = topologies[uid]
        if topology not in passed_uids_grouped.keys(): #If not in new dict - add
            passed_uids_grouped[topology] = [uid]
        else:
            passed_uids_grouped[topology].append(uid) #Otherwise append

    with open(logfile, 'a+') as file: #append to file and create it if it does not exist
        file.write('PDB filter\n')
        file.write('#uids before:'+str(len(sequences))+'\t')
        file.write('#uids after:'+str(len(passed_uids))+'\t')
        file.write('%uids passed:'+str(np.around(len(passed_uids)/len(sequences)*100, decimals = 1))+'\n')
        file.write('#Topologies before:'+str(num_topologies_before)+'\t')
        file.write('#Topologies after:'+str(len(passed_uids_grouped))+'\t')
        file.write('%Topologies passed:'+str(np.around(len(passed_uids_grouped)/num_topologies_before*100, decimals = 1))+'\n')

	
	
    return passed_uids_grouped

def overx(passed_uids_grouped, outdir, sequences, H_groups):
	'''Write all uids to their respective topology in fasta format
	'''

	#Write all passed groups to newline separated
	#text files in partitions of 1000
	topology_groupings = open('topology_groupings.txt', 'w')

	x=2 #At least two entries per topology for a pair
	#write passed_uids_grouped_over_x to fasta
	for group in passed_uids_grouped:
		group_dir = outdir+'fasta/'+group
		uids = passed_uids_grouped[group]
		if len(uids)<2: #Need to uids to make a pair
			continue
		#shuffle uids to make sure the groupings are random
		random.Random(2).shuffle(uids)
		used_hgroups = [] #Save one uid per H-group
		selected_uids = [] #The chosen uids
		for uid in uids:
			hgroup = H_groups[uid]
			if hgroup in used_hgroups:#Only use each h-group once
				continue
			else:
				used_hgroups.append(hgroup)
				selected_uids.append(uid)
		if len(selected_uids)<2: #If there is not at least one pair
			continue
		os.mkdir(group_dir)
		topology_groupings.write(group_dir+'\n')
		if len(selected_uids)%2 != 0: #If even pairs can not be created. Exclude one.
			selected_uids = selected_uids[1:]
		for uid in selected_uids:
			sequence = sequences[uid]
			with open(group_dir+'/'+uid+'.fa', "w") as file:
				file.write('>'+uid+'\n')
				j = 0 #index
				while j<len(sequence):
					file.write(sequence[j:j+60]+'\n')
					j+=60
	topology_groupings.close()

	return None

#####MAIN#####
args = parser.parse_args()

H_groups, topologies = read_tsv(args.H_groups[0]) #Read H-groups and uids
sequences = read_fasta(args.sequences[0]) #Read the clustered sequences
failed_pdb_filter = pd.read_csv(args.failed_pdb_filter[0], sep = '\n', header = None) #Get uids that failed the pdb filter
outdir = args.outdir[0]
logfile = args.logfile[0]

#Number of topologies before
grouped_seqs_by_topology = get_groups(topologies, sequences)
#Assess which passed
passed_uids_grouped = get_passed_uids(failed_pdb_filter, sequences, topologies, len(grouped_seqs_by_topology)) #Get uids that passed the pdb filter and group them
overx(passed_uids_grouped, outdir, sequences, H_groups)
