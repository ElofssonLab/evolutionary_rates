#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import numpy as np
import pandas as pd
import subprocess
from Bio import pairwise2
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Parse dssp and match to both structural and sequence alignments''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to input directory. include / in end''')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory. include / in end''')
parser.add_argument('fastadir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to fasta directory. include / in end''')
parser.add_argument('group', nargs=1, type= str,
                  default=sys.stdin, help = '''group.''')
parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')

#Max acc surface areas for each amino acid according to empirical measurements in:
#Tien, Matthew Z et al. “Maximum allowed solvent accessibilites of residues in proteins.”
#PloS one vol. 8,11 e80635. 21 Nov. 2013, doi:10.1371/journal.pone.0080635
max_acc = { 'A':121,
			'R':265,
			'N':187,
			'D':187,
			'C':148,
			'E':214,
			'Q':214,
			'G':97,
			'H':216,
			'I':195,
			'L':191,
			'K':230,
			'M':203,
			'F':228,
			'P':154,
			'S':143,
			'T':163,
			'W':264,
			'Y':255,
			'V':165,
			'X':192 #Average of all other maximum surface accessibilites
		  }



#FUNCTIONS
def read_fasta(filepath):
	'''Read fasta sequence
	'''


	with open(filepath, 'r') as file:
		sequence = ''
		for line in file:
			line = line.rstrip() #remove \n
			if line[0] == '>':
				uid = line[1:]
			else:
				sequence += line

	return(sequence)

def run_dssp(indir,group, uid):
	'''Run dssp
	'''
	DSSP = '/home/pbryant/dssp'
	command = DSSP +' '+indir+group+'/'+uid+'.pdb'
	outp = subprocess.check_output(command, shell = True)#run dssp
	outp = outp.decode()
	with open(indir+group+'/dssp/'+uid+'.pdb.dssp', 'w') as file:
		file.write(outp)
	return None

def match_dssp_to_aln(df, indir, outdir, fastadir):
	'''Match alignments to dssp surface acc and 2ndary str
	'''
	DSSP = '/home/p/pbryant/pfs/dssp'
	dssp_dict = {}
	fasta_dict = {}

	#Create dict to save results
	results_dict = {'ss1_seqaln':[],'acc1_seqaln':[],
	                'ss2_seqaln':[], 'acc2_seqaln':[],
	                'ss1_straln':[], 'acc1_straln':[],
                    'ss2_straln':[],'acc2_straln':[]
                    }


	for index, row in df.iterrows():

		group = row['group']
		uid1 = row['uid1']
		uid2 = row['uid2']

        #Get dssp and fasta info
		if uid1 not in dssp_dict.keys():
			dssp_file = indir+group+'/dssp/'+uid1+'.pdb.dssp'
			if not os.path.isfile(dssp_file):
				run_dssp(indir,group, uid1)

			sequence1, secondary_str1, surface_acc1 = parse_dssp(indir+group+'/dssp/'+uid1+'.pdb.dssp')
			dssp_dict[uid1] = [sequence1, secondary_str1, surface_acc1]
			fasta_dict[uid1] = read_fasta(fastadir+group+'/'+uid1+'.fa')
		if uid2 not in dssp_dict.keys():
			dssp_file = indir+group+'/dssp/'+uid2+'.pdb.dssp'
			if not os.path.isfile(dssp_file):
                                run_dssp(indir,group, uid2)
			sequence2, secondary_str2, surface_acc2 = parse_dssp(indir+group+'/dssp/'+uid2+'.pdb.dssp')
			dssp_dict[uid2] = [sequence2, secondary_str2, surface_acc2]
			fasta_dict[uid2] = read_fasta(fastadir+group+'/'+uid2+'.fa')

		#For sequence and structure alignments
		for suffix in ['_seqaln', '_straln']:

			#Match sequences to alignments
			aln1 = row['seq1'+suffix]
			aln2 = row['seq2'+suffix]
			#Save gapless alignments
			gapless_aln1 = ''
			gapless_aln2 = ''

			for i in range(len(aln1)):
				if aln1[i] == '-' or aln2[i] == '-':
					    continue
				else:
					gapless_aln1 += aln1[i]
					gapless_aln2 += aln2[i]

			(ss1, acc1) = match(gapless_aln1, dssp_dict[uid1], fasta_dict[uid1])
			(ss2, acc2) = match(gapless_aln2, dssp_dict[uid2], fasta_dict[uid2])

			results_dict['ss1'+suffix].append(ss1)
			results_dict['acc1'+suffix].append(acc1)

			results_dict['ss2'+suffix].append(ss2)
			results_dict['acc2'+suffix].append(acc2)

	#Save to df
	for key in results_dict:
		df[key] = results_dict[key]

	#Write new df to outdir
	df.to_csv(outdir+group+'_df.csv')
	return None

def parse_dssp(filepath):
	'''Parse dssp output and return sequence, secondary structure and surface accessibilities.
	'''

	secondary_str = [] #save str
	surface_acc = [] #save acc
	sequence = '' #save residues
	fetch_lines = False #Don't fetch unwanted lines
	try:
		with open(filepath, 'r') as file:
			for line in file:
				if fetch_lines == True:
					line = line.rstrip()
					residue = line[13]
					str_i = line[16]
					acc_i = line[35:38].strip()

					if residue == '!' or residue == '!*':
						continue
					else:
					#Normalize acc_i by the max acc surface area for the specific amino acid
					#Round to whole percent
						acc_i_norm = round((float(acc_i)/max_acc[residue])*100, )
						acc_i_norm = min(acc_i_norm, 100) #Should not be over 100 percent

						#Save
						secondary_str.append(str_i)
						surface_acc.append(acc_i_norm)
						sequence = sequence + residue
				if '#' in line:
					fetch_lines = True
					#now the subsequent lines will be fetched
	except:
		if os.path.isfile(filepath):
			raise ValueError('File', filepath,'exists but could not be read')
		else:
			raise ValueError('No file', filepath)
	return(sequence, secondary_str, surface_acc)

def match(gapless_aln, info, org_seq):
	'''Extract secondary structure annotations and surface acc matching the aligned sequences'''

	dssp_seq = info[0]
	all_ss = info[1]
	all_acc = info[2]

	matched_ss = ''
	matched_acc = ''

	#align to original sequence
	aln1 = pairwise2.align.globalxx(org_seq, gapless_aln)
	aln2 = pairwise2.align.globalxx(org_seq, dssp_seq)

	if len(aln1)<1 or len(aln2)<1:
		raise IOError('Alignments of length 0')

	seq1 = aln1[0][1] #aln to org
	seq2 = aln2[0][1] #dssp to org
	dsspi=0
	for i in range(len(seq1)):
		if seq1[i] != '-':
			if seq2[i] != '-':
				matched_ss += all_ss[dsspi]+','
				matched_acc += str(all_acc[dsspi])+','
			else:
				matched_ss += '-'+','
				matched_acc += '-'+','
		if seq2[i] != '-':
			dsspi += 1

	return (matched_ss, matched_acc)

#MAIN
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
fastadir = args.fastadir[0]
group = args.group[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)
df = df[df['group']==group]
if len(df) ==0:
	raise IOError('Empty')
match_dssp_to_aln(df, indir, outdir, fastadir)
