#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import pdb
import os

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that filters entries in CATH
        according to a filter file containing PDB ids. All entries in each topology are then 
	randomly paired and written to individual folders. One folder per topology.
	Each toplogoical folder has subfolders, each with at most 100 pairs.''')

parser.add_argument('--H_groups', nargs=1, type= str,
                  default=sys.stdin, help = 'path to H_groups.')
parser.add_argument('--sequences', nargs=1, type= str,
                  default=sys.stdin, help = 'path to clustered sequences.')
parser.add_argument('--failed_pdb_filter', nargs=1, type= str,
                  default=sys.stdin, help = 'path to uids that failed the pdb filter.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'path output directory.')



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


def get_passed_uids(failed_pdb_filter, sequences, topologies):
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
    passed_uids_grouped = {} #Sequences grouped by H-group
    for uid in passed_uids:
        topology = topologies[uid]
        if topology not in passed_uids_grouped.keys(): #If not in new dict - add
            passed_uids_grouped[topology] = [uid]
        else:
            passed_uids_grouped[topology].append(uid) #Otherwise append


    return passed_uids_grouped

def overx(passed_uids_grouped, outdir, sequences):
    '''Write all uids to their respective topology in fasta format
    '''

    #Write all passed groups to newline separated
    #text files in partitions of 1000
    topology_groupings = open('topology_groupings.txt', 'w')

    x=2 #At least two entries per topology for a pair
    #write passed_uids_grouped_over_x to fasta
    for group in passed_uids_grouped:
        group_dir = outdir+'/fasta/'+group
        uids = passed_uids_grouped[group]
        if len(uids)<2:
            continue
        os.mkdir(group_dir)
        for i in range(0,len(uids),200):
            uids_200 = uids[i:i+200]
            subgroup_dir = group_dir+'/'+str(i)
            os.mkdir(subgroup_dir)
            topology_groupings.write(subgroup_dir+'\n')		
            for uid in uids_200:
                sequence = sequences[uid]
                with open(subgroup_dir+'/'+uid+'.fa', "w") as file:
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

grouped_sequences = get_groups(topologies, sequences) #Group sequences
passed_uids_grouped = get_passed_uids(failed_pdb_filter, sequences, topologies) #Get uids that passed the pdb filter and group them
overx(passed_uids_grouped, outdir, sequences)
