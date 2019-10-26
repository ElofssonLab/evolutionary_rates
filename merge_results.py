#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import glob
import subprocess
from collections import Counter
import os
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that collects results and merges them
                                                into a unified datframe. A selection on how much of
                                                the shortest sequence in each pair is made according
                                                to a threshold.''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to input directory.')
parser.add_argument('threshold', nargs=1, type= float,
                  default=sys.stdin, help = 'Threshold for percent aligned of the shortest sequence in each pair.')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'path output directory.')
parser.add_argument('dssp_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to dssp.''')
parser.add_argument('fastadir', nargs=1, type= str,
                default=sys.stdin, help = '''path to fastadir.''')
parser.add_argument('gitdir', nargs=1, type= str,
                default=sys.stdin, help = '''path to gitdir.''')


def tsv_to_df(globpath, dfname):
    '''Get results from tsv files into
    '''
    #Get results from structure alignments
    all_files = glob.glob(globpath)     # advisable to use os.path.join as this makes concatenation OS independent
    df_from_each_file = [pd.read_csv(f, sep='\t') for f in all_files]
    for dataframe, filename in zip(df_from_each_file, all_files):
    	hgroup = filename.split('/')[-3]
    	dataframe['H_group'] = hgroup
    	dataframe['C.'] = hgroup.split('.')[0]+'.'
    	dataframe['C.A.'] = hgroup.split('.')[0]+'.'+hgroup.split('.')[1]
    	df = pd.concat(df_from_each_file, ignore_index=True)
    	df.to_csv(dfname)
    return df

def get_sequence_alignments(globpath, dfname):
    '''Get all sequence alignments
    '''

    #Save all alignements and their metrics for the hhalign sequence alignments
    all_files = glob.glob(globpath)

    uid1 = []
    uid2 = []
    l1 = []
    s1 = []
    e1 = []
    l2 = []
    s2 = []
    e2 = []
    aln_len = []
    identities = []
    e_values = []
    probabilities = []
    seq1 = []
    seq2 = []
    groups = []
    for name in all_files:
        groups.append(name.split('/')[-2])
        with open(name, 'r') as file:
            for line in file:
                line = line.rstrip() #remove \n
                if '>' in line:
                    uid = line[1:8]
                    line = line.split('|') #Split
                    if len(uid1)>len(uid2): #append to uid lists
                        uid2.append(uid)
                        l2.append(int(line[1].split('=')[1]))
                        s2.append(int(line[2].split('=')[1]))
                        e2.append(int(line[3].split('=')[1]))
                    else:
                        uid1.append(uid)
                        l1.append(int(line[1].split('=')[1]))
                        s1.append(int(line[2].split('=')[1]))
                        e1.append(int(line[3].split('=')[1]))
                        aln_len.append(int(line[4].split('=')[1]))
                        identities.append(float(line[5].split('=')[1]))
                        e_values.append(float(line[6].split('=')[1]))
                        probabilities.append(float(line[7].split('=')[1]))
                else:
                    sequence = line
                    if len(seq1)>len(seq2):
                        seq2.append(sequence)
                    else:
                        seq1.append(sequence)


    seq_aln_df = pd.DataFrame(list(zip(uid1, uid2, aln_len, identities, e_values, probabilities, seq1, seq2, l1, l2, s1, s2, e1, e2)), columns = ['uid1','uid2','aln_len','identity', 'e_value', 'probability', 'seq1', 'seq2', 'l1', 'l2', 's1', 's2', 'e1', 'e2'])
    seq_aln_df.to_csv(dfname)

    return seq_aln_df


def get_structure_alignments(globpath, dfname):
    '''Get all sequence alignments
    '''

    #Save all alignements and their metrics for the TMalign sequence alignments
    all_files = glob.glob(globpath)
    uid1 = []
    uid2 = []
    l1 = []
    l2 = []
    aln_len = []
    identities = []
    seq1 = []
    seq2 = []
    groups = []
    for name in all_files:
        groups.append(name.split('/')[-3])
        with open(name, 'r') as file:
            for line in file:
                line = line.rstrip() #remove \n
                if '>' in line:
                    uid = line[1:8]
                    line = line.split('|') #Split
                    if len(uid1)>len(uid2): #append to uid lists
                        uid2.append(uid)
                        l2.append(int(line[1].split('=')[1]))

                    else:
                        uid1.append(uid)
                        l1.append(int(line[1].split('=')[1]))
                        aln_len.append(int(line[2].split('=')[1]))
                        identities.append(float(line[3].split('=')[1]))
                else:
                    sequence = line
                    if len(seq1)>len(seq2):
                        seq2.append(sequence)
                    else:
                        seq1.append(sequence)


    str_aln_df = pd.DataFrame(list(zip(uid1, uid2, aln_len, identities,seq1, seq2, l1, l2)), columns = ['uid1','uid2','aln_len','identity', 'seq1', 'seq2', 'l1', 'l2'])
    str_aln_df.to_csv(dfname)

    return str_aln_df

def get_lddt(globpath, dfname):
    '''Get all lddt results into a df
    '''
    all_files = glob.glob(globpath)

    uid1 = []
    uid2 = []
    lddt_scores = []
    for name in all_files:
        with open(name, 'r') as file:
            uids = name.split('/')[-1].split('.')[0].split('_')
            uid1.append(uids[0])
            uid2.append(uids[1])
            for row in file:
                if 'Global LDDT score:' in row:
                    row = row.split()
                    score = float(row[-1])
                    lddt_scores.append(score)
                    break

    lddt_df = pd.DataFrame(list(zip(uid1, uid2, lddt_scores)), columns = ['uid1','uid2','lddt_scores'])
    lddt_df.to_csv(dfname)

    return lddt_df

def percent_aligned(df):
    '''Calculate percent of shortest sequence in each pair that has been aligned
    '''

    percent_aligned = []
    for index, row in df.iterrows():
        min_len = min(row['l1'], row['l2'])
        percent_aligned.append(row['aln_len']/min_len)

    df['percent_aligned'] = percent_aligned

    return df

def create_df(results_path, t, dssp, fastadir, gitdir):
    '''Gets results and parses them into a unified dataframe
    '''

    #Get results from sequence alignments into df
    globpath = results_path+'*/TMscore/*seq.tsv'
    seq_df_name = results_path+'/from_seq_df.csv'
    from_seq_df = tsv_to_df(globpath, seq_df_name)

    #Get results from structure alignments into df
    globpath = results_path+'*/TMalign/*str.tsv'
    str_df_name = results_path+'/from_str_df.csv'
    from_str_df = tsv_to_df(globpath, str_df_name)

    #get sequence alignments
    globpath = results_path+'*/*.aln'
    seq_aln_name = results_path+'/seq_aln_df.csv'
    seq_aln_df = get_sequence_alignments(globpath, seq_aln_name)

    #get sequence alignments from structure alignments
    globpath = results_path+'*/TMalign/*.aln'
    str_aln_name = results_path+'/str_aln_df.csv'
    str_aln_df = get_structure_alignments(globpath, str_aln_name)

    #Get lddt results for sequence alignments
    globpath = results_path+'*/TMscore/*.lddt'
    seq_lddt_name = results_path+'/seq_lddt_df.csv'
    seq_lddt_df = get_lddt(globpath, seq_aln_name)

    #Get lddt results for structure alignments
    globpath = results_path+'*/TMalign/*.lddt'
    str_lddt_name = results_path+'/str_lddt_df.csv'
    str_lddt_df = get_lddt(globpath, str_aln_name)

    #Merge seq aln dfs
    complete_seq_df = pd.merge(from_seq_df, seq_aln_df, on=['uid1', 'uid2'], how='left')
    complete_seq_df = pd.merge(seq_lddt_df, complete_seq_df, on=['uid1', 'uid2'], how='left')
    complete_seq_df = complete_seq_df.dropna() #Drop NANs
    complete_seq_df = percent_aligned(complete_seq_df) #Calculate percent aligned
    complete_seq_df.to_csv(results_path+'/complete_seq_df.csv')

    #Merge str aln dfs
    complete_str_df = pd.merge(from_str_df, str_aln_df, on=['uid1', 'uid2'], how='left')
    complete_str_df = pd.merge(str_lddt_df, complete_str_df, on=['uid1', 'uid2'], how='left')
    complete_str_df = complete_str_df.dropna() #Drop NANs
    complete_str_df = percent_aligned(complete_str_df) #Calculate percent aligned
    complete_str_df.to_csv(results_path+'/complete_str_df.csv')

    #Select entries from complete seq df based on percent aligned
    complete_seq_df_t = complete_seq_df[complete_seq_df['percent_aligned']>=t]

    #Merge sequence and str dataframes after percent aligned selection
    complete_df = pd.merge(complete_seq_df_t, complete_str_df, on=['uid1', 'uid2'], how='left')
    complete_df = complete_df.dropna() #Drop NANs

    #Rename x with seqaln and y with straln
    cols = ['lddt_scores', 'MLAAdist', 'RMSD', 'aln_len', 'identity', 'seq1', 'seq2', 'l1', 'l2', 'percent_aligned', ]
    for col in cols:
        complete_df = complete_df.rename(columns={col+'_x': col+'_seqaln', col+'_y': col+'_straln'})

    #rename H_group_x to H_group
    complete_df = complete_df.rename(columns={'H_group_x': 'H_group'})
    #Save df
    complete_df.to_csv(results_path+'/complete_df.csv')
    #Run DSSP - better done in parallel
    hgroups = [*Counter(complete_df['H_group']).keys()]

    os.mkdir(results_path+'/dssp_dfs/') #Make dssp dir
    for group in hgroups:
        command = gitdir+'match_dssp.py '+results_path+' '+results_path+'/dssp_dfs/ '+fastadir+' '+group+' '+results_path+'/complete_df.csv'
        outp = subprocess.check_output(command, shell = True)#run dssp

    #Consolidate dssp dfs
    all_files = glob.glob(results_path+'/dssp_dfs/*.csv')
    df_from_each_file = [pd.read_csv(f) for f in all_files]
    complete_dssp_df = pd.concat(df_from_each_file, ignore_index=True)
    complete_dssp_df = complete_dssp_df.dropna() #Drop NANs
    complete_dssp_df.to_csv(results_path+'/complete_df.csv')


    #Calculate contacts - better done in parallel
    hgroups = [*Counter(complete_dssp_df['H_group']).keys()]
    os.mkdir(results_path+'/contacts/') #Make contact dir

    for group in hgroups:
        command = gitdir+'contact_calculations.py '+results_path+' '+results_path+'/contacts/ '+fastadir+' '+results_path+'/complete_dssp_df.csv'+' '+group
        outp = subprocess.check_output(command, shell = True)#run dssp

    #Consolidate contact dfs
    all_files = glob.glob(results_path+'/contacts/*.csv')
    df_from_each_file = [pd.read_csv(f) for f in all_files]
    complete_dssp_contact_df = pd.concat(df_from_each_file, ignore_index=True)
    complete_dssp_contact_df = complete_dssp_contact_df.dropna() #Drop NANs
    complete_dssp_contact_df.to_csv(results_path+'/complete_df.csv')


    #Reduce cardinalities
    return None

#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
t = args.threshold[0]
outdir = args.outdir[0]
dssp = args.dssp_path[0]
fastadir = args.fastadir[0]
gitdir = args.gitdir[0]
create_df(indir, t, dssp, fastadir, gitdir)
