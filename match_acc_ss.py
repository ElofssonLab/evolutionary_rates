#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pdb

def match_acc(df, suffix):

    DIFF_ACC = [] #Save fractional ids

    s1 = ''
    s2 = ''
    for index, row in df.iterrows():

        acc1 = row['acc1'+suffix].split(',')[:-1]
        acc2 = row['acc2'+suffix].split(',')[:-1]
        diff = 0 #Difference btw sequences in terms of buried or exposed residues

        for i in range(len(acc1)):
            if acc1[i] == '-' or acc2[i] == '-': #If gap - continue
                continue
            else:
                if int(acc1[i]) > 25:
                    s1 = 'E'
                else:
                    s1 = 'B'

                if int(acc2[i]) > 25:
                    s2 = 'E'
                else:
                    s2 = 'B'

                if s1 == s2: #both buried/exposed
                    continue
                else:
                    diff +=1

        if len(acc1) == 0:
            pdb.set_trace()
        DIFF_ACC.append(diff/len(acc1))

    return DIFF_ACC

def match_ss(df, suffix):

    IDSS = [] #Save fractional ids
    helics = ['G', 'H', 'I']
    strands = ['E', 'S']
    loops = ['T', 'S', 'B']

    ss = {'G':'H',
          'H':'H',
          'I':'H',
          'E':'S',
          'B':'S',
          'T':'L',
          'S':'L',
          ' ':'L' #COILS
         }


    for index, row in df.iterrows():
         #States - reset
        states = {
            'HH':0,
            'SS':0,
            'LL':0,
            'HS':0,
            'SH':0,
            'HL':0,
            'LH':0,
            'SL':0,
            'LS':0
            }

        ss1 = row['ss1'+suffix]
        ss1 = ss1.split(',')[:-1]
        ss2 = row['ss2'+suffix]
        ss2 = ss2.split(',')[:-1]

        for i in range(len(ss1)):
            if ss1[i] == '-' or ss2[i] == '-':
                continue
            else:
                match = ss[ss1[i]]+ss[ss2[i]]
                states[match] += 1

        #Calculate IDSS
        T = (states['HH']+states['SS']+states['LL'])
        N = T + states['HS']+states['SH']+states['HL']+states['LH']+states['SL']+states['LS']
        if N != 0:
            IDSS.append(T/N)
        else:
            pdb.set_trace()

    return IDSS
