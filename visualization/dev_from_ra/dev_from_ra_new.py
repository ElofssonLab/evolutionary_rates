#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse
from scipy import stats
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols
import time
'''
When using statsmodels in scientific publication, please consider using the following citation:
Seabold, Skipper, and Josef Perktold. “Statsmodels: Econometric and statistical modeling with python.”
Proceedings of the 9th Python in Science Conference. 2010.
'''

#Custom
from percent_ss import parse_ss
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that calculates deviations for running averages.''')

parser.add_argument('--topdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--hgroupdf', nargs=1, type= str,
default=sys.stdin, help = 'path to df.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either median or average.')

parser.add_argument('--avdf', nargs=1, type= str,
default=sys.stdin, help = 'Dataframe with averages.')

def dev_from_av(avdf, df, score, aln_type, cardinality, max_seqdist):
    '''Calculate avearge dev from total running average within group and significance
    '''

    if cardinality == '_AA20':
            cardinality = ''

    absolute_deviations = [] #absolute of all devs
    rco_extreme = []
    avs = [] #Save average score
    avdevs = [] #Save deviations
    step = 0.1


    mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
    start = np.round(min(mldists),2)
    start = float(str(start)[0:3])
    end = min(max(mldists), max_seqdist) #End at max_seqdist
    scores = np.asarray(df[score+aln_type])
    js = [] #Save js

    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_df = df[df['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]

        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        if len(below_df)<1: #May be discontinuous in step
            continue
        cut_scores = np.asarray(below_df[score+aln_type])
        if calc == 'average':
            av= np.average(cut_scores)
        if calc == 'median':
            av= np.median(cut_scores)
        avs.append(av)

        x = np.round(j-step/2, 2)
        js.append(x)
        tav = avdf[avdf['ML  distance']==x][score+aln_type].values[0] #total average in current interval
        avdevs.append(av-tav)

        #Plot deviation against RCO
        absolute_deviations.extend(np.absolute(cut_scores-tav))
        rco1 = [*below_df['RCO1']]
        rco2 = [*below_df['RCO2']]


    #Do a t-test to calculate if the deviation is significant for each H-group
    #The variances will be unequal, as each H-group only encompasses very feq points
    #True values = 0, no deviation from total average
    truevals = np.zeros(len(df)) #size should not matter since zeros
    statistic, pvalue = stats.ttest_ind(avdevs, truevals, equal_var = False)
    #plt.plot(js, avs)

    return np.average(avdevs), pvalue, js, avs

def plot_partial(partial_df, partial_merged, avdf, name, score, aln_type, cardinality, title, results_dir, color):
    '''RA plots of partial dfs and the total RA
    '''
    topologies = [*partial_df['Topology']]
    ylims = {'RMSD':[0,4], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6], 'lddt_scores': [0.2,1.0], 'DIFFC':[0,1], 'TMscore': [0.2,1.0]}
    grad_ylims = {'RMSD':[-0.1,0.1], 'lddt_scores':[-0.025, 0.025], 'DIFFSS':[-0.025, 0.025], 'DIFF_ACC':[-0.025, 0.025], 'DIFFC':[-0.04, 0.04], 'TMscore':[-0.025, 0.025]}
    mldists = [*partial_df[score+aln_type+'_seqdists']]
    scores = [*partial_df[score+aln_type+'_ra']]
    stds = np.array(avdf[score+aln_type+'_std']) #std dev for total ra
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    matplotlib.rcParams.update({'font.size': 7})
    #Plot RA per topology
    #Percent within std
    percent_within_std = pd.DataFrame(list(zip(np.arange(0.05,6,0.1), np.zeros(60),  np.zeros(60))),columns = ['seqdist', 'within','total'])
    percent_within_std['seqdist'] = np.round(percent_within_std['seqdist'], 2)
    for i in range(len(partial_df)):
        top = topologies[i]
        ax.plot(mldists[i],scores[i], alpha = 0.1, color = color, linewidth =1)
        #Get percentage within 1 std
        for j in range(len(mldists[i])):
            sel = avdf[avdf['ML  distance']==mldists[i][j]]
            s = sel[score+aln_type].values[0]
            std = sel[score+aln_type+'_std'].values[0]
            d = scores[i][j]
            ind = percent_within_std[percent_within_std['seqdist']==mldists[i][j]].index[0]
            percent_within_std['total'][ind]+=1 #Add count
            if d > s-std and d <s+std: #If within std
                percent_within_std['within'][ind]+=1 #Add count


    #Plot total RA for topologies
    step = 0.1
    start = 0
    end = 6
    total_top_js = []
    total_top_ra = []
    cardinality = ''
    for j in np.arange(start+step,end+step,step):
        if np.round(j, 2) == end: #Make sure to get endpoints
            below_df = partial_merged[partial_merged['MLAAdist'+cardinality+aln_type]<=j]
        else:
            below_df = partial_merged[partial_merged['MLAAdist'+cardinality+aln_type]<j]

        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
        if len(below_df)<1: #May be discontinuous in step
            continue
        cut_scores = np.asarray(below_df[score+aln_type])

        total_top_ra.append(np.average(cut_scores))
        total_top_js.append(np.round(j-step/2,2))

    ax.plot(total_top_js, total_top_ra, color = color, linewidth = 1.5, label = 'Topology', alpha = 0.7)
    ax.plot(avdf['ML  distance'], avdf[score+aln_type], color = 'g', linewidth = 1.5, label = 'Broad dataset')
    #plot stddev
    ax.plot(avdf['ML  distance'], np.array(avdf[score+aln_type])+np.array(stds), '--', c = 'g', linewidth = 1) #positive stds
    ax.plot(avdf['ML  distance'], np.array(avdf[score+aln_type])-np.array(stds), '--', c = 'g', linewidth = 1, label = 'Standard deviation') #negative stds
    ax.legend()
    ax.set_title(title)
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_ylim(ylims[score])
    ax.set_xlabel('ML AA20 distance')
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score')
    else:
        ax.set_ylabel(score)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(results_dir+name, format = 'png')

    # #Scatterplot
    # fig = plt.figure(figsize=(11,11)) #set figsize
    # plt.scatter(partial_merged['MLAAdist'+cardinality+aln_type],partial_merged[score+aln_type],alpha = 0.2, color = 'b', s = 1)
    # plt.plot(total_top_js, total_top_ra, color = 'b', linewidth = 3, label = 'Topology')
    # plt.plot(avdf['ML  distance'], avdf[score+aln_type], color = 'r', linewidth = 3, label = 'Broad dataset')
    #
    # plt.legend()
    # plt.xlim([0,9.1])
    # plt.xticks([0,1,2,3,4,5,6,7,8,9])
    # plt.ylim(ylims[score])
    # plt.xlabel('ML AA20 distance')
    # if score == 'lddt_scores':
    #     plt.ylabel('lDDT score')
    # else:
    #     plt.ylabel(score)
    # fig.savefig(results_dir+'scatter_'+name, format = 'png')

    #Plot gradients
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    #T-test
    partial_avdf = avdf[avdf['ML  distance']<6]
    statistic, pvalue = stats.ttest_ind(np.gradient(total_top_ra), np.gradient(partial_avdf[score+aln_type]), equal_var = False)
    ax.plot(total_top_js, np.gradient(total_top_ra), color = 'b', linewidth = 1.5, label = 'Topology\npval:'+str(np.round(pvalue,2)))
    ax.plot(avdf['ML  distance'], np.gradient(avdf[score+aln_type]), color = 'r', linewidth = 1.5, label = 'Broad dataset')
    ax.legend()
    ax.set_xlim([0,9.1])
    ax.set_xticks([0,1,2,3,4,5,6,7,8,9])
    ax.set_ylim(grad_ylims[score])
    ax.set_xlabel('ML AA20 distance')
    if score == 'lddt_scores':
        ax.set_ylabel('lDDT score gradients')
    else:
        ax.set_ylabel(score+' gradients')
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(results_dir+'gradients_'+name, format = 'png')

    #Close plots to avoid to many being open at the same time
    plt.close()

    return percent_within_std

def class_percentages(df):
    '''Calculate class percentages
    '''

    percentages = []
    for C in [1.,2.,3.,4.]:
        try:
            percentages.append(len(df[df['C._x']==C])/len(df))
        except:
            percentages.append(0)

    return percentages

def anova(top_metrics_sig , features, aln_type, results_dir, colors):
    '''Perform anova and plot violinplots
    '''
    #Add size
    features.append(score+aln_type+'_sizes')

    #Make violinplots
    titles = {'RCO':'RCO', 'aln_len'+aln_type:'Aligned length', 'l':'Length', 'percent_aligned'+aln_type:'% Aligned',
    'K':'KR','D':'DE','Y':'YWFH','T':'TSQN','C':'CVMLIA', 'P':'PG', '-':'Gap', 'L':'Loop', 'S':'Sheet', 'H':'Helix',
    score+aln_type+'classes':'Class', score+aln_type+'_sizes':'Group size'}
    ylabels = {'RCO':'Average RCO', 'aln_len'+aln_type:'Aligned length', 'l':'Length', 'percent_aligned'+aln_type:'% Aligned',
    'K':'%','D':'%','Y':'%','T':'%','C':'%', 'P':'%', '-':'%', 'L':'%', 'S':'%', 'H':'%'}
    matplotlib.rcParams.update({'font.size': 7})
    for feature in features:
        fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
        sns.violinplot(data = top_metrics_sig, x = 'Significance', y = feature+'_av', palette=colors, order=['+','Non', '-'])
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel('Topology Average')
        ax.set_title(titles[feature])
        fig.tight_layout()
        fig.savefig(results_dir+'volin_'+feature+'_'+score+aln_type+'.png', format = 'png')

        plt.close()

    #Plot class distributions
    features.append(score+aln_type+'_classes')
    class_freq = pd.DataFrame()
    cs = []
    sigs = []
    percs = []
    for sig in ['+', 'Non', '-']:
        sel = top_metrics_sig[top_metrics_sig['Significance']==sig]
        counts = Counter(sel[score+aln_type+'_classes_av'])
        cs.extend(['α', 'β', 'αβ', 'Few SS'])
        sigs.extend([sig]*4)
        for key in [1.,2.,3.,4.]:
            try:
                percs.append(counts[key]/len(sel))
            except:
                percs.append(0)

    class_freq['Significance'] = sigs
    class_freq['Percentage'] = percs
    class_freq['Class'] = cs
    fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
    sns.barplot(data=class_freq, x='Class',y='Percentage',hue='Significance', palette = colors)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Percentage')
    ax.set_title('Class')
    fig.tight_layout()
    fig.savefig(results_dir+'volin_class_'+score+aln_type+'.png', format = 'png')



    #Calculate fraction retained
    eta_sq = []
    omega_sq = []
    Pr = []


    #Fit ANOVA, using all features simultaneously
    anova_cols = []
    #Rename gap cols
    top_metrics_sig = top_metrics_sig.rename(columns={'-_av':'gap_av'})
    for feature in features:
        if feature == '-':
            feature = 'gap'
        test_str = feature+'_av ~ C(Significance)'
        results = ols(test_str, data=top_metrics_sig).fit()
        aov_table = sm.stats.anova_lm(results, typ=2)
        outp = anova_table(aov_table)
        eta_sq.append(outp.iloc[0]['eta_sq'])
        omega_sq.append(outp.iloc[0]['omega_sq'])
        Pr.append(outp.iloc[0]['PR(>F)'])
        anova_cols.append(feature)

    anova_df =  pd.DataFrame(data = {'Feature': anova_cols, 'eta_sq': eta_sq, 'omega_sq':omega_sq, 'PR(>F)':Pr})
    anova_df.to_csv(results_dir+'anova_tests.csv')
    return anova_df

def anova_table(aov):
    aov['mean_sq'] = aov[:]['sum_sq']/aov[:]['df']

    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])

    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*aov['mean_sq'][-1]))/(sum(aov['sum_sq'])+aov['mean_sq'][-1])

    cols = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[cols]
    return aov

def ttest_features(df, catdf, score, aln_type):
    '''Perform t-tests for features - looking into differences btw groups
    '''

    features = ['RCO', 'aln_len'+aln_type, 'l', 'percent_aligned'+aln_type]

    results = {}
    for feature in features:
        if feature == 'RCO':
            x = np.concatenate((np.array(df['RCO1']), np.array(df['RCO2'])))
            y = np.concatenate((np.array(catdf['RCO1']), np.array(catdf['RCO2'])))
        elif feature == 'l':
            x = np.concatenate((np.array(df['l1'+aln_type]), np.array(df['l2'+aln_type]))) #the lengths for the sequence alignments represents the picked up consensus from hhsearch
            y = np.concatenate((np.array(catdf['l1'+aln_type]), np.array(catdf['l2'+aln_type])))
        else:
            x = df[feature]
            y = catdf[feature]

        statistic, pvalue = stats.ttest_ind(x,y, equal_var = False)
        #Z-scores
        z = (np.average(x)-np.average(y))/(np.std(x)/np.sqrt(len(df)))
        #Average
        av = np.average(x)
        results[feature] = [statistic, pvalue, z, av]

    #Compare AA6 distributions and secondary structure elements on topology group level
    for key in ['P', 'C', 'K', 'T', 'D', 'Y', '-', 'L', 'S', 'H' ]:
        x = np.concatenate((np.array(df[key+'1'+aln_type]), np.array(df[key+'2'+aln_type]))) #the lengths for the sequence alignments represents the picked up consensus from hhsearch
        y = np.concatenate((np.array(catdf[key+'1'+aln_type]), np.array(catdf[key+'2'+aln_type])))
        statistic, pvalue = stats.ttest_ind(x,y, equal_var = False)
        #Z-scores
        z = (np.average(x)-np.average(y))/(np.std(x)/np.sqrt(len(df)))
        #Average
        av = np.average(x)
        results[key] = [statistic, pvalue, z, av]
    return results

def percent_sig_in_set(pos_sig, nonsig_df, neg_sig, features, results_dir, aln_type, perc_keys, colors):
    '''Calculate % sig in each set for different features
    Divide into pos, neg and non-deviating
    '''
    matplotlib.rcParams.update({'font.size': 7})
    total = len(pos_sig)+ len(nonsig_df) + len(neg_sig)
    x = np.arange(3)
    w = 1/4 #width of bar
    titles = {'RCO':'RCO', 'aln_len'+aln_type:'Aligned length', 'l':'Length', 'percent_aligned'+aln_type:'% Aligned',
    'K':'KR','D':'DE','Y':'YWFH','T':'TSQN','C':'CVMLIA', 'P':'PG', '-':'-'}
    for key in features:
        #Pos
        pos = pos_sig[pos_sig[key+'_pval']<0.05/total]
        pos_pos = len(pos[pos[key+'_tstat']>0])
        pos_neg = len(pos[pos[key+'_tstat']<0])
        pos_non = len(pos_sig)-(pos_pos+pos_neg)
        #Neg
        neg = neg_sig[neg_sig[key+'_pval']<0.05/total]
        neg_pos = len(neg[neg[key+'_tstat']>0])
        neg_neg = len(neg[neg[key+'_tstat']<0])
        neg_non = len(neg_sig)-(neg_pos+neg_neg)
        #Non
        non = nonsig_df[nonsig_df[key+'_pval']<0.05/total]
        non_pos = len(non[non[key+'_tstat']>0])
        non_neg = len(non[non[key+'_tstat']<0])
        non_non = len(nonsig_df)-(non_pos+non_neg)

        fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
        ax.bar(x, [100*pos_pos/len(pos_sig), 100*pos_non/len(pos_sig), 100*pos_neg/len(pos_sig)],w, label='+', color = colors[0])
        ax.bar(x+w, [100*non_pos/len(nonsig_df), 100*non_non/len(nonsig_df), 100*non_neg/len(nonsig_df)],w, label='Non', color = colors[1])
        ax.bar(x+2*w, [100*neg_pos/len(neg_sig), 100*neg_non/len(neg_sig), 100*neg_neg/len(neg_sig)],w, label='-', color = colors[2])
        ax.set_title(titles[key])
        ax.set_xticks(x+w)
        ax.set_xticklabels(['Pos', 'Non', 'Neg'])
        ax.set_xlabel('Feature significance')
        ax.set_ylabel('% in each set')
        ax.set_ylim([0,100])
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.tight_layout()
        plt.legend()
        fig.savefig(results_dir+key+'_'+score+aln_type+'.png', format = 'png')
        plt.close()

    #Plot size distributions
    fig, ax = plt.subplots(figsize=(4.5/2.54,4.5/2.54))
    sns.distplot(pos_sig[score+aln_type+'_sizes'],  label='+', hist = False, color = colors[0])
    sns.distplot(nonsig_df[score+aln_type+'_sizes'],  label='Non', hist = False, color = colors[1])
    sns.distplot(neg_sig[score+aln_type+'_sizes'], label='-', hist = False, color = colors[2])
    plt.legend()
    plt.xscale("log")
    ax.set_ylim([0,0.02])
    ax.set_yticks(np.arange(0,0.025,0.005))
    ax.set_xticks([10,100,1000])
    ax.set_xlabel('log Group size')
    ax.set_ylabel('Density')
    ax.set_title('Group size')
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(results_dir+'sizes_'+score+aln_type+'.png', format = 'png')
    plt.close()


    #Look at H, L, S

    #plt.show()
    return None

def AA6_distribution(df, aln_type):
    '''Calculate AA6 distribution in sequence (gaps have their own state as well)
    Groupings: K=[KR],D=[DE],Y=[YWFH],T=[TSQN],C=[CVMLIA], P=[PG], "-"="-"
    '''
    #Save percentages:
    percentage_dict = {
    'K1'+aln_type:[],
    'K2'+aln_type:[],
    'D1'+aln_type:[],
    'D2'+aln_type:[],
    'Y1'+aln_type:[],
    'Y2'+aln_type:[],
    'T1'+aln_type:[],
    'T2'+aln_type:[],
    'C1'+aln_type:[],
    'C2'+aln_type:[],
    'P1'+aln_type:[],
    'P2'+aln_type:[],
    '-1'+aln_type:[],
    '-2'+aln_type:[]
    }


    #reset index
    df = df.reset_index()
    t1 = time.time()
    for i in range(len(df)):
        row = df.iloc[i]
        l1 = len(row['AA6_1_straln'])
        l2 = len(row['AA6_2_straln'])
        s1 = Counter(row['AA6_1_straln'])
        s2 = Counter(row['AA6_2_straln'])

        #Get s1 ad s2 percentages
        for key in ['P', 'C', 'K', 'T', 'D', 'Y', '-']:
            try:
                percentage_dict[key+'1'+aln_type].append(s1[key]/l1)
            except:
                percentage_dict[key+'1'+aln_type].append(0)
            try:
                percentage_dict[key+'2'+aln_type].append(s2[key]/l2)
            except:
                percentage_dict[key+'2'+aln_type].append(0)

    print('Calculating AA6',time.time()-t1,'s')

    for key in percentage_dict:
        df[key] = percentage_dict[key]

    return df, percentage_dict.keys()

def percent_ss(df):
    '''Calculate % H,L and S in aligned structures
    '''

    return None


def three_sets_comparison(catdf_s, top_metrics, score, aln_type, cardinality, features,  perc_keys, outdir):
    '''Compares positively, negatively and non-deviating groups in their
    deviation from the total running average.
    '''
    colors = ['royalblue', 'cadetblue', 'cornflowerblue']
    try:
        os.mkdir(outdir+score+aln_type)
    except:
        print('Dir '+outdir+score+aln_type+'/'+' exists')
    #Plot size against average deviation
    matplotlib.rcParams.update({'font.size': 7})
    fig, ax = plt.subplots(figsize=(12/2.54,6/2.54))
    ax.scatter(top_metrics['lddt_scores_straln_sizes_av'],top_metrics['lddt_scores_straln_av_dev'],s=0.5)
    ax.set_ylim([-0.2,0.2])
    ax.set_xticks(np.arange(0,5000,2000))
    ax.set_xlabel('Group size')
    ax.set_ylabel('Deviation')
    ax.set_title('Size and Deviation')
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+score+aln_type+'/'+'size_dev_'+score+aln_type+'.png', format = 'png')
    plt.close()

    #Get significant
    sig_df = top_metrics[top_metrics[score+aln_type+'_ra_pval']<0.05/len(top_metrics)]
    #Get pos deviating>0
    pos_sig = sig_df[sig_df[score+aln_type+'_av_dev']>0]
    pos_sig['Significance'] = '+'
    #Get neg deviating<0
    neg_sig = sig_df[sig_df[score+aln_type+'_av_dev']<0]
    neg_sig['Significance'] = '-'
    #check
    if len(pos_sig)+len(neg_sig) != len(sig_df):
        pdb.set_trace()

    #Get non significant
    nonsig_df = top_metrics[top_metrics[score+aln_type+'_ra_pval']>=0.05/len(top_metrics)]
    nonsig_df['Significance'] = 'Non'
    #Get data from cat_df matching sig topologies
    pos_sig_merged = pd.merge(pos_sig, catdf_s, left_on='Topology', right_on='group', how='left')
    neg_sig_merged = pd.merge(neg_sig, catdf_s, left_on='Topology', right_on='group', how='left')
    #Get data from cat_df matching sig topologies
    nonsig_df_merged = pd.merge(nonsig_df, catdf_s, left_on='Topology', right_on='group', how='left')

    #Calculate ANOVA
    top_metrics_sig = pd.concat([pos_sig, neg_sig])
    top_metrics_sig = pd.concat([top_metrics_sig, nonsig_df])
    anova_df = anova(top_metrics_sig, features, aln_type, outdir+'/'+score+aln_type+'/', colors)

    #Calculate percentages of sig for each feature in each set
    #percent_sig_in_set(pos_sig, nonsig_df, neg_sig, features, outdir+score+aln_type+'/', aln_type,  perc_keys, colors)


    #Plot all RAs per top group
    top_metrics_merged = pd.merge(top_metrics, catdf_s, left_on='Topology', right_on='group', how='left')
    total_within_std = plot_partial(top_metrics,top_metrics_merged, avdf, score+aln_type+'_ra_per_top.png', score, aln_type, cardinality, 'All Topologies with 10', outdir+'/'+score+aln_type+'/', 'b')

    #Plot the RAs of the pos and neg sig groups
    pos_within_std = plot_partial(pos_sig,pos_sig_merged, avdf, score+aln_type+'_ra_pos_sig.png', score, aln_type, cardinality, 'Pos. set', outdir+'/'+score+aln_type+'/', colors[0])
    non_within_std = plot_partial(nonsig_df, nonsig_df_merged, avdf, score+aln_type+'_ra_non_sig.png', score, aln_type, cardinality, 'Non. set', outdir+'/'+score+aln_type+'/', colors[1])
    neg_within_std = plot_partial(neg_sig, neg_sig_merged, avdf, score+aln_type+'_ra_neg_sig.png', score, aln_type, cardinality, 'Neg. set', outdir+'/'+score+aln_type+'/', colors[2])

    #Plot percent within std
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    ax.plot(total_within_std['seqdist'],100*total_within_std['within']/total_within_std['total'], label = 'Total', color = 'g', linewidth = 2)
    #ax.plot(pos_within_std['seqdist'],pos_within_std['within']/pos_within_std['total'], label = 'Set A', color = colors[0], linewidth = 2)
    #ax.plot(non_within_std['seqdist'],non_within_std['within']/non_within_std['total'], label = 'Set B', color = colors[2], linewidth = 2)
    #ax.plot(neg_within_std['seqdist'],neg_within_std['within']/neg_within_std['total'], label = 'Set C', color = colors[1], linewidth = 2)
    #ax.set_xlabel('ML AA20 Distance')
    ax.set_ylabel('%')
    ax.set_title('Within 1 Std')
    ax.set_xticks(np.arange(0,9.1,1))
    ax.set_xlabel('ML AA20 distance')
    ax.set_ylim([60,100])
    ax.set_yticks([60,70,80,90,100])
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+score+aln_type+'/'+'within_std'+score+aln_type+'.png', format = 'png')
    plt.close()

    #Concat
    cat_dev = pd.concat([pos_sig_merged, neg_sig_merged])
    classes = ['Mainly Alpha', 'Mainly Beta', 'Alpha Beta', 'Few SS']
    print(score+aln_type)
    #print('%pos sig', class_percentages(pos_sig_merged))
    #print('%neg sig', class_percentages(neg_sig_merged))
    #print('%nonsig', class_percentages(nonsig_df_merged))
    cat_dev = pd.concat([cat_dev, nonsig_df_merged])



    return None

#####MAIN#####
args = parser.parse_args()
topdf = pd.read_csv(args.topdf[0])
hgroupdf = pd.read_csv(args.hgroupdf[0])
outdir = args.outdir[0]
calc = args.calc[0]
avdf = pd.read_csv(args.avdf[0])

cardinality = '_AA20'

#Get topology from hgroupdf
tops = []
hgroups = [*hgroupdf['group']]
for hg in hgroups:
    hg = hg.split('.')
    tops.append(hg[0]+'.'+hg[1]+'.'+hg[2])

hgroupdf['C.A.T.'] = tops
#rename col
hgroupdf = hgroupdf.rename(columns={'group':'H_group'})
hgroupdf = hgroupdf.rename(columns={'C.A.T.':'group'})

#rename TMscore cols
hgroupdf = hgroupdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})
topdf = topdf.rename(columns={'TMscore':'TMscore_seqaln', 'TMscore_high':'TMscore_straln'})

catdf = pd.concat([topdf, hgroupdf])
#Rename class column
catdf = catdf.rename(columns={'C._x':'Class'})

#The ones should actually be zeros
catdf['RCO1']=catdf['RCO1'].replace([1], 0)
catdf['RCO2']=catdf['RCO2'].replace([1], 0)
topcounts = Counter(catdf['group'])
vals = np.array([*topcounts.values()])
num_tops_with10 = len(np.where(vals>9)[0]) #Get all topologies with at least 10 values
print('Fraction of topologies with at least 10 entries: '+str(num_tops_with10)+'/'+str(len(vals)),num_tops_with10/len(vals))
topologies = np.array([*topcounts.keys()])
topologies = topologies[np.where(vals>9)[0]]

#Save pvalues
top_metrics = pd.DataFrame()
top_metrics['Topology'] = topologies
cardinality = '' #AA20

for score in ['lddt_scores', 'TMscore', 'DIFFC', 'RMSD', 'DIFFSS', 'DIFF_ACC']:
    for aln_type in ['_straln', '_seqaln']:


        #select below 6 using seq or str
        catdf_s = catdf[catdf['MLAAdist'+aln_type]<=6]
        if score == 'lddt_scores':
            #Save statistical assessment - this only has to be done for 1 score (but both alntypes) - will be the same in groups
            stat_results = {}
            features = ['RCO', 'aln_len'+aln_type, 'l', 'percent_aligned'+aln_type,'P', 'C', 'K', 'T', 'D', 'Y',
            '-', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix
            for key in features:
                stat_results[key] = np.zeros((555,4))
            catdf_s, perc_keys = AA6_distribution(catdf_s, aln_type) #Get AA6 frequencies
            catdf_s = parse_ss(catdf_s, aln_type) #Get % ss
        avs_from_line = [] #save avs from line and pvals
        pvals = []
        all_js = []
        all_avs = []
        gradients = []
        toplens = []
        av_RCO = []
        sizes = [] #group sizes
        classes = []

        for top in topologies:
            df = catdf_s[catdf_s['group']==top]
            toplens.append(len(df))
            av_from_line, pvalue, js, avs = dev_from_av(avdf, df, score, aln_type, cardinality, 6)
            avs_from_line.append(av_from_line)
            pvals.append(pvalue)
            all_js.append(js)
            all_avs.append(avs)
            gradients.append(np.gradient(avs))
            sizes.append(len(df))
            classes.append(df['Class'].values[0]) #append class

            if score == 'lddt_scores':
                #Test for deviating features within group
                results = ttest_features(df, catdf_s, score, aln_type)
                for key in results: #each key is a feature
                    stat_results[key][len(all_avs)-1] = results[key] #statisic, pvalue, z


        if score == 'lddt_scores':
            for key in features:
                top_metrics[key+'_tstat'] = stat_results[key][:,0]
                top_metrics[key+'_pval'] = stat_results[key][:,1]
                top_metrics[key+'_z'] = stat_results[key][:,2]
                top_metrics[key+'_av'] = stat_results[key][:,3]
        top_metrics[score+aln_type+'_ra_pval'] = pvals
        top_metrics[score+aln_type+'_av_dev'] = avs_from_line
        top_metrics[score+aln_type+'_seqdists'] = all_js
        top_metrics[score+aln_type+'_ra'] = all_avs
        top_metrics[score+aln_type+'_gradients'] = gradients
        top_metrics[score+aln_type+'_sizes_av'] = sizes
        top_metrics[score+aln_type+'_classes_av'] = classes

        #plt.scatter(top_metrics['lddt_scores_straln_sizes'], top_metrics['lddt_scores_straln_av_dev'], s= 5)
        #sel = top_metrics[top_metrics['lddt_scores_straln_sizes']<500]
        #plt.scatter(sel['lddt_scores_straln_sizes'], sel['lddt_scores_straln_av_dev'], s= 5)
        #Make plots
        three_sets_comparison(catdf_s, top_metrics, score, aln_type, cardinality, features, perc_keys, outdir)





#top_metrics.to_csv(outdir+'top_metrics.csv')
