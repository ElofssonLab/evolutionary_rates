#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.spatial import distance
from Bio import pairwise2
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate the relative contact order for pdb files.''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to input directory. include / in end''')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory. include / in end''')
parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')





#FUNCTIONS

def calculate_rco(df, indir, outdir):
	'''Calculate the relative contact order
	'''


	contact_dict = {}
	fasta_dict = {}

	groups = [*Counter(df[*['group']]).keys()] #Get unique groups
	uids = [*Counter(uid1+uid2).keys()] #Get unique uids


	for group in groups:
		uid1 = df[df['group']==group]['uid1']
		uid2 = df[df['group']==group]['uid2']
		uids = [*Counter(uid1+uid2).keys()] #Get unique uids


		for uid in uids:
			contacts, sequence = read_cbs(indir+group+'/'+uid+'.pdb')
			contact_dict[uid1] = [contacts, sequence]



			C = 0 #Keep track of the number of conserved contacts in the alignment
			for i in range(len(mc1)):
				c1 = mc1[i]
				c2 = mc2[i]
				for j in c1:
					if j in c2:#If the contacts at the same position is shared.
						C+=1
			try:
				diff = 1-(C/(M+N-C))
			except:
				diff = 'NaN'
			if suffix == '_seqaln':
				DIFFC_seqaln.append(diff)
			else:
				DIFFC_straln.append(diff)

	#Set new columns in df
	df['DIFFC_seqaln'] = DIFFC_seqaln
	df['DIFFC_straln'] = DIFFC_straln
	#Write new df to outdir
	df.to_csv(outdir+group+'_df.csv')
	return None

def read_cbs(pdbfile):
	'''Get the C-betas from a pdb file.
	'''
	three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
	sequence = ''
	pos = [] #Save positions in space
	prev_res = -1 #Keep track of potential alternative residues
	with open(pdbfile, 'r') as file:
		for line in file:
			record = parse_atm_record(line)
			if record['atm_name'] == 'CB':
				if record['res_no'] == prev_res:
					continue
				else:
					prev_res = record['res_no']
					pos.append(np.array([record['x'], record['y'], record['z']]))
					sequence += three_to_one[record['res_name']]
			if record['atm_name'] == 'CA' and record['res_name'] == 'GLY':
				prev_res = record['res_no']
				pos.append(np.array([record['x'], record['y'], record['z']]))
				sequence += three_to_one[record['res_name']]
	contacts = get_contacts(pos)
	return contacts, sequence


def parse_atm_record(line):

	record = defaultdict()
	record['name'] = line[0:6].strip()
	record['atm_no'] = int(line[6:11])
	record['atm_name'] = line[12:16].strip()
	record['res_name'] = line[17:20].strip()
	record['chain'] = line[21]
	record['res_no'] = int(line[22:26])
	record['insert'] = line[26].strip()
	record['resid'] = line[22:29]
	record['x'] = float(line[30:38])
	record['y'] = float(line[38:46])
	record['z'] = float(line[46:54])
	record['occ'] = float(line[54:60])
	record['B'] = float(line[60:66])

	return record

def get_contacts(pos):
	contacts = [] #Save each residue's contacts
	for i in range(len(pos)):
		contacts.append([])
		for j in range(i+5, len(pos)):
			dist = distance.euclidean(pos[i], pos[j])
			if dist < 8:
				contacts[i].append(j)


	return contacts

#MAIN
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)

contact_order(df, indir, outdir)
