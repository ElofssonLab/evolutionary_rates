#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

from sklearn.model_selection import train_test_split
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
                  default=sys.stdin, help = 'If optimization is to be performed or not (Bool).')

#FUNCTIONS

def create_features(df):
    '''Get features
    '''

    double_features = ['CD', 'RCO', 'l', 'P', 'C', 'K', 'T', 'D', 'Y', '-', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix, contact density

    #Get features
    all_features = ['MLAAdist_straln', 'percent_aligned_straln']
    X = []
    X.append(np.array(df['MLAAdist_straln']))
    X.append(np.array(df['percent_aligned_straln']))
    for x in double_features:
        if x == 'RCO' or x == 'CD':
            for i in ['1', '2']:
                X.append(np.array(df[x+i]))
                all_features.append(x+i)
        else:
            for i in ['1', '2']:
                X.append(np.array(df[x+i+'_straln']))
                all_features.append(x+i+'_straln')
    X = np.array(X) #Convert to np array
    X=X.T #Transpose

    #Get deviation
    y = np.asarray(df['lddt_scores_straln_dev']) #binned_rmsds

    return(X, y, all_features)

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


def plot_predictions(X, y, all_features):
    '''Predict and plot
    '''
    matplotlib.rcParams.update({'font.size': 7})
    #Make train and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    #Fit
    #'classify__min_samples_split': 2, 'classify__n_estimators': 500, 'classify__verbose': 0} -0.033 +/- 0.11
    rfreg = RandomForestRegressor(n_estimators=500,
                        	    verbose=0,
                                max_depth= None,
                                min_samples_split=2)

    rfreg.fit(X_train,y_train)

    #Feature importance
    importances = rfreg.feature_importances_
    imp_df = pd.DataFrame()
    imp_df['Feature'] = all_features
    imp_df['Importance'] = [*importances]
    fig, ax = plt.subplots(figsize=(12/2.54,12/2.54))
    sns.barplot(y="Feature", x="Importance", data=imp_df)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir+'feature_importances.png', format = 'png')
    plt.close()


    #predict
    rfreg_predictions = rfreg.predict(X_test)
    #Average error
    error = rfreg_predictions-y_test
    print(np.average(np.absolute(error)))
    #Plot ED against error
    make_kde(X_test[:,0], error, 'AA20 ED', 'Error (lDDT score)', [0,6],[-0.1,0.1], outdir+'ed_vs_error.png', False)
    make_kde(rfreg_predictions,y_test, 'Pred. dev. (lDDT score)', 'True dev. (lDDT score)',[-0.1,0.1], [-0.1,0.1], outdir+'pred_vs_true.png', True)


def make_kde(x,y, xlabel, ylabel, xlim, ylim, outname, get_R):
    '''Makes a kdeplot and saves it
    '''
    fig, ax = plt.subplots(figsize=(12/2.54,12/2.54))
    sns.kdeplot(x,y, shade = True, kind = 'kde', cmap = 'Blues')
    if get_r == True:
        print(outname,pearsonr(x,y)[0])
        ax.plot(xlim,ylim, 'darkblue')
        plt.annotate(str(pearsonr(x,y)[0]), (-0.5,0.5))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

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

#Assign data and labels
X, y, all_features = create_features(df)

#Grid search and cross validate
param_grid = {'classify__n_estimators': [750,1000],
        'classify__max_depth': [None],
        'classify__min_samples_split': [2]}


#RandomForestRegressor
rfreg = RandomForestRegressor()
pipe = Pipeline(steps=[('classify', rfreg)])
if optimize == True:
    pdb.set_trace()
    best_clf = parameter_optimization(param_grid, pipe, X, y)
else:
    plot_predictions(X, y, all_features)
