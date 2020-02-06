#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A random RandomForestClassifier for predicting
                                                lddt between structural alignments.''')

parser.add_argument('--dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')



#FUNCTIONS

def create_features(df):
    '''Get features
    '''

    double_features = ['CD', 'RCO', 'l', 'P', 'C', 'K', 'T', 'D', 'Y', '-', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix, contact density

    #Get features
    X = []
    X.append(np.array(df['MLAAdist_straln']))
    X.append(np.array(df['percent_aligned_straln']))
    for x in double_features:
        if x == 'RCO' or x == 'CD':
            for i in ['1', '2']:
                X.append(np.array(df[x+i]))
        else:
            for i in ['1', '2']:
                X.append(np.array(df[x+i+'_straln']))
    X = np.array(X) #Convert to np array
    X=X.T #Transpose



    #Get deviation
    y = np.asarray(df['lddt_scores_straln_dev']) #binned_rmsds

    return(X, y)

def make_kde(x,y, xlabel, ylabel, outname):
    '''Makes a kdeplot and saves it
    '''
    fig, ax = plt.subplots(figsize=(12/2.54,12/2.54))
    sns.kdeplot(x,y, shade = True, kind = 'kde', cmap = 'Blues')
    print(pearsonr(x,y))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.savefig(outname, format = 'png')
    plt.close()

#MAIN
args = parser.parse_args()
df = pd.read_csv(args.dataframe[0])
outdir = args.outdir[0]
#Assign data and labels
X,y = create_features(df)

#Make train and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
#RandomForestClassifier
rfreg = RandomForestRegressor(n_estimators=100, bootstrap = True, max_features = 'sqrt')
# Fit on training data
rfreg.fit(X_train, y_train)
#predict
rfreg_predictions = rfreg.predict(X_test)
#Average error
error = rfreg_predictions-y_test
print(np.average(np.absolute(error)))
#Plot ED against error
make_kde(X_test[:,0], error, 'ED', 'Error (lDDT score)', outdir+'ed_vs_error.png')
make_kde(rfreg_predictions,y_test, 'Pred. dev.', 'True dev.', outdir+'pred_vs_true.png')
