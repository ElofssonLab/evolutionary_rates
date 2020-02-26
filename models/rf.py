#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A random RandomForestClassifier for predicting
                                                lddt between structural alignments.''')

parser.add_argument('--dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

parser.add_argument('--optimize', nargs=1, type= int,
                  default=sys.stdin, help = 'If optimization (1) is to be performed or not (0) (Bool).')

parser.add_argument('--get_all', nargs=1, type= int,
                  default=sys.stdin, help = 'If all features (1) or only sequence features (0) are to be used (Bool).')

#FUNCTIONS

def create_features(df, get_all):
    '''Get features
    '''

    nicer_names = {'CD':'Contact density', 'RCO':'RCO', 'l':'Length', 'P':'PG', 'C':'CVMLIA',
    'K':'KR','D':'DE','Y':'YWFH','T':'TSQN', '-':'Gap', 'L':'Loop', 'S':'Sheet', 'H':'Helix'}

    sequence_features = ['l', 'P', 'C', 'K', 'T', 'D', 'Y', '-']
    structure_features = ['CD', 'RCO', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix, contact density
    if get_all == True:
        double_features = sequence_features+structure_features
    else:
        double_features = sequence_features
    #Get features
    all_features = ['AA20 ED', '% Aligned']
    X = []
    X.append(np.array(df['MLAAdist_straln']))
    X.append(np.array(df['percent_aligned_straln']))
    for x in double_features:
        if x == 'RCO' or x == 'CD':
            for i in ['1', '2']:
                X.append(np.array(df[x+i]))
                all_features.append(nicer_names[x]+i)
        else:
            for i in ['1', '2']:
                X.append(np.array(df[x+i+'_straln']))
                all_features.append(nicer_names[x]+' '+i)
    X = np.array(X) #Convert to np array
    X=X.T #Transpose

    #Get deviation
    y = np.asarray(df['lddt_scores_straln_dev']) #binned_rmsds

    #Get true lDDT scores
    z_true = np.asarray(df['lddt_scores_straln'])
    #Average lDDT scores
    z_av = np.asarray(df['lddt_scores_straln_av'])

    return(X, y, z_true, z_av, all_features)

def parameter_optimization(param_grid, pipe, X, y):
    '''Optimize parameters
    '''
    #Scores to optimize for
    #Scores to optimize for
    scores = ['neg_mean_absolute_error']
    best_classifiers = {}
    store_means = {}
    store_stds = {}

    #Optimize parameters
    for score in scores:
        with open(outdir+'opt.txt', 'a+') as file:
            file.write("# Tuning hyper-parameters \n")

            #Grid search and cross validate
            #cv determines the cross-validation splitting strategy, 5 specifies the number of folds (iterations) in a (Stratified)KFold
            clf = GridSearchCV(pipe,
                               param_grid= param_grid, cv=5, #cv = 5, Stratified Kfold
                               scoring='neg_mean_absolute_error',
                               n_jobs=-1)#use all cores
            #Fit on train_data and write to file
            clf.fit(X,y)
            #Write to file
            file.write("Best parameters set found on development set:" + '\n')
            file.write(str(clf.best_params_))
            file.write('\n' + '\n' + "Grid scores on development set:" + '\n')
            means = clf.cv_results_['mean_test_score']
            store_means[score] = means
            stds = clf.cv_results_['std_test_score']
            store_stds[score] = stds
            for mean, std, params in zip(means, stds, clf.cv_results_['params']):
                file.write('mean test score: ')
                file.write("%0.3f (+/-%0.03f) for %r"
                      % ( mean, std * 2, params) + '\n')

        return clf.best_estimator_


def plot_predictions(X, y, z_true, z_av, all_features):
    '''Predict and plot
    '''
    matplotlib.rcParams.update({'font.size': 7})
    #Store metrics
    error = []
    error_av = []
    R_dev = []
    R_lDDT = []
    R_av = []
    if os.path.exists('pred0.csv'): #If predictions have been made
        pred_df = pd.read_csv('pred0.csv')
        imp_df = pd.read_csv('feature_imp.csv')
    else:
        #No previous predictions
        #Clf
        #'classify__min_samples_split': 2, 'classify__n_estimators': 500, 'classify__verbose': 0} -0.033 +/- 0.11
        rfreg = RandomForestRegressor(n_estimators=500,
                                    verbose=5,
                                    max_depth= None,
                                    min_samples_split=2,
                                    n_jobs=-1)
        #Make train and test set
        cv = KFold(n_splits=5, random_state=42, shuffle=True)
        #Save feature importances
        importances = []
        imp_df = pd.DataFrame()
        imp_df['Feature'] = all_features

        for i, (train_index, test_index) in enumerate(cv.split(X)):
            #Fit
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            rfreg.fit(X_train,y_train)
            #Feature importance
            importances.append(rfreg.feature_importances_)
            #predict
            rfreg_predictions = rfreg.predict(X_test)
            #Prediction df
            pred_df = pd.DataFrame()
            pred_df['pred'] = rfreg_predictions
            pred_df['true'] = y_test
            pred_df['ED']= X_test[:,0]
            pred_df['lDDT_true'] = z_true[test_index]
            pred_df['lDDT_av'] = z_av[test_index]
            pred_df.to_csv('pred'+str(i)+'.csv')
            #Error in deviation
            error.append(np.average(np.absolute(rfreg_predictions-y_test)))
            R_dev.append(pearsonr(rfreg_predictions,y_test)[0])
            #Error in lDDT
            R_lDDT.append(pearsonr(z_true[test_index],z_av[test_index]+rfreg_predictions)[0])
            #Error to average
            R_av.append(pearsonr(z_true[test_index],z_av[test_index])[0])
            error_av.append(np.average(np.absolute(z_true[test_index]-z_av[test_index])))

        importances = np.array(importances) #Convert to array
        av_imp = []
        std_imp = []
        for i in range(len(all_features)):
            av_imp.append(np.average(importances[:,i]))
            std_imp.append(np.std(importances[:,i]))

        #Calculate average
        imp_df['Importance'] = av_imp
        imp_df['std'] = std_imp
        imp_df = imp_df.sort_values(['Importance'], ascending = False).reset_index(drop=True)
        #Save
        imp_df.to_csv('feature_imp.csv')

    #Average and std
    if len(error)==0: #If already predicted
        for i in range(5):
            df = pd.read_csv('pred'+str(i)+'.csv')
            rfreg_predictions = df['pred']
            y_test = df['true']
            lDDT_true = df['lDDT_true']
            lDDT_av = df['lDDT_av']
            error.append(np.average(np.absolute(rfreg_predictions-y_test)))
            R_dev.append(pearsonr(rfreg_predictions,y_test)[0])
            R_lDDT.append(pearsonr(lDDT_true,lDDT_av+rfreg_predictions)[0])

    print('Error of model:',np.round(np.average(error),3),'+/-', np.round(np.std(error),5))
    #print('Error to RA:',np.round(np.average(error_av),3),'+/-', np.round(np.std(error_av),5))
    print('R_dev:',np.round(np.average(R_dev),3),'+/-', np.round(np.std(R_dev),3))
    print('R_lDDT:',np.round(np.average(R_lDDT),3),'+/-', np.round(np.std(R_lDDT),3))
    #print('R_av:',np.round(np.average(R_av),3),'+/-', np.round(np.std(R_av),3))

    #Plot
    fig, ax = plt.subplots(figsize=(12/2.54,12/2.54))
    sns.barplot(y="Feature", x="Importance", data=imp_df)
    plt.errorbar(imp_df["Importance"],np.arange(0,len(all_features)), yerr=None, xerr=imp_df['std'], fmt='.')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'feature_importances.png', format = 'png')
    plt.close()



    #Plot ED against error
    make_kde(pred_df['ED'], np.array(pred_df['pred']-pred_df['true']), 'AA20 ED', 'Error (lDDT score)', [0,6],[-0.1,0.1], np.arange(0,6), np.arange(-0.1,0.11, 0.05), outdir+'ed_vs_error.png', False, 0, (0,0))
    make_kde(pred_df['pred'],pred_df['true'], 'Pred. dev. (lDDT score)', 'True dev. (lDDT score)',[-0.1,0.1], [-0.1,0.1], np.arange(-0.1,0.11, 0.05), np.arange(-0.1,0.11, 0.05), outdir+'pred_vs_true_dev.png', True, 0.72, (-0.08,0.05))
    make_kde(pred_df['lDDT_av']+pred_df['pred'], pred_df['lDDT_true'], 'Pred. lDDT score', 'True lDDT score',[0.2,1], [0.2,1], np.arange(0.2,1.1, 0.1), np.arange(0.2,1.1, 0.1), outdir+'pred_vs_true_lddt.png', True, 0.92, (0.25,0.7))
    make_kde(pred_df['lDDT_av'], pred_df['lDDT_true'], 'Average lDDT score', 'True lDDT score',[0.2,1], [0.2,1], np.arange(0.2,1.1, 0.1), np.arange(0.2,1.1, 0.1), outdir+'av_vs_true_lddt.png', True, 0.84, (0.25,0.7))


def make_kde(x,y, xlabel, ylabel, xlim, ylim, xticks, yticks, outname, get_R, R, cords):
    '''Makes a kdeplot and saves it
    '''
    fig, ax = plt.subplots(figsize=(6/2.54,6/2.54))
    sns.kdeplot(x,y, shade = True, shade_lowest = False, kind = 'kde', cmap = 'Blues')
    if get_R == True:
        print(outname,pearsonr(x,y)[0])
        ax.plot(xlim,ylim, 'darkblue')
        plt.annotate('Pearson R: '+str(R), cords)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_xticks(xticks)
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outname, format = 'png')
    plt.close()

#MAIN
args = parser.parse_args()
df = pd.read_csv(args.dataframe[0])
outdir = args.outdir[0]
optimize = bool(args.optimize[0])
get_all = bool(args.get_all[0])

#Assign data and labels
X, y, z_true, z_av, all_features = create_features(df, get_all)

#Grid search and cross validate
param_grid = {'classify__n_estimators': [750,1000],
        'classify__max_depth': [None],
        'classify__min_samples_split': [2]}

#RandomForestRegressor
rfreg = RandomForestRegressor()
pipe = Pipeline(steps=[('classify', rfreg)])
if optimize == True:
    best_clf = parameter_optimization(param_grid, pipe, X, y)
else:
    plot_predictions(X, y, z_true, z_av, all_features)
