#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Extract information from .eps figure.''')

parser.add_argument('figure', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to train dataframe in .csv.')



def parse_eps(fig):
        '''Extract points from eps
        '''

        data_dict = {}
        i = 0
        get_points = False
        with open(fig, 'r') as figfile:
            for line in figfile:
                line = line.rstrip()
                if line == "/bg { 1 0 0 rgb } def" or line == '/bg { 0 1 0 rgb } def':
                    i+=1
                    data_dict[i] = [[],[]]
                    get_points = True
                    continue
                if get_points == True:
                    line = line.split()
                    if len(line)<6:
                        continue
                    if line[4] == 'cl':
                            get_points = False
                            continue
                    if line[2] == '1.79':
                        data_dict[i][0].append(float(line[0]))
                        data_dict[i][1].append(float(line[1]))




        #Need to rescale to go from 0.1-6 and 0.2-5
        #ML dist
        x1 = np.array(data_dict[3][0])
        x2 = np.array(data_dict[4][0])
        x = np.concatenate([np.array(data_dict[3][0]), np.array(data_dict[4][0])])
        print(min(x),max(x))
        x_new = (x-min(x))
        x1_new = x1-min(x)
        x1_new = x1_new*6/max(x_new)
        x2_new = x2-min(x)
        x2_new = x2_new*6/max(x_new)
        #x_new = x_new+0.1
        print(min(x_new),max(x_new))

        #RMSD: 0.2-4.2
        y1= np.array(data_dict[3][1])
        y2 = np.array(data_dict[4][1])
        y = np.concatenate([np.array(data_dict[3][1]), np.array(data_dict[4][1])])
        print(min(y),max(y))
        y_new = (y-min(y))
        y1_new = y1-min(y)
        y1_new = y1_new*5.9/max(y_new)
        y1_new = y1_new+0.2
        y2_new = y2-min(y)
        y2_new = y2_new*5.9/max(y_new)
        y2_new = y2_new+0.2

        print(min(y_new),max(y_new))

        #Plot
        fig, ax = plt.subplots(figsize=(8/2.54,8/2.54)) #set figsize
        matplotlib.rcParams.update({'font.size': 7})
        ax.scatter(x1_new, y1_new, s= 0.6, marker = 's', c = 'r')
        ax.scatter(x2_new, y2_new, s= 0.6, marker = 's', c = 'lime')
        ax.set_ylim([-0.1,4.1])
        ax.set_yticks([0,1,2,3,4])
        ax.set_xticks([0,1,2,3,4,5,6])
        ax.set_xlabel('ML AA20 distance')
        ax.set_ylabel('RMSD')
        fig.tight_layout()
        fig.savefig('extracted.png', format = 'png')
        x = np.concatenate([x1_new, x2_new])
        y = np.concatenate([y1_new, y2_new])

        return x,y



#MAIN
args = parser.parse_args()
fig = args.figure[0]

x,y = parse_eps(fig)
data = pd.DataFrame()
data['seqdist']=x
data['RMSD']=y
data.to_csv('eps_df.csv')
