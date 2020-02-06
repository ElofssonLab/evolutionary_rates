#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from collections import Counter

from sklearn.linear_model import LinearRegression

from model_inputs import split_on_h_group
import matplotlib.pyplot as plt
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A random RandomForestClassifier for predicting
                                                lddt between structural alignments.''')

parser.add_argument('--dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('--out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')



#FUNCTIONS

def create_features(df):
    '''Get features
    '''

    single_features = ['percent_aligned'+aln_type]
    double_features = ['CD', 'RCO', 'l', 'P', 'C', 'K', 'T', 'D', 'Y', '-', 'L', 'S', 'H'] #L,S,H = loop, sheet, helix, contact density

    #Get features
    X = []
    X.append('MLAAdist_straln')
    X.append(np.array('percent_aligned_straln'))
    for x in double_features:
        if x == 'RCO' or x == 'CD':
            for i in ['1', '2']:
                X.append(np.array(df[x+i]))
        else:
            for i in ['1', '2']:
                X.append(np.array(df[x+i+'_straln']))
    X = np.array(X) #Convert to np array


    #Get deviation
    y = np.asarray(df['lddt_scores_straln_dev']) #binned_rmsds

    return(X, y)

#MAIN
args = parser.parse_args()
df = pd.read_csv(args.dataframe[0])
out_dir = args.out_dir[0]
#Assign data and labels
X,y = create_features(df)


#RandomForestClassifier
rfreg = RandomForestRegressor(n_estimators=100, bootstrap = True, max_features = 'sqrt')
# Fit on training data
clf.fit(X_train, y_train)
rfreg.fit(X_train, y_train)
#predict
clf_predictions = clf.predict(X_valid)
rfreg_predictions = rfreg.predict(X_valid)
#Average error
average_error = np.average(np.absolute(clf_predictions-y_valid))
print(average_error)
average_error = np.average(np.absolute(rfreg_predictions-y_valid))
print(average_error)
