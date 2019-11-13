#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import svgutils.transform as sg
from svgutils.compose import *
import sys
import argparse

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that combines svg figures into publication ready svgs.''')

parser.add_argument('--fig1', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 1.')

parser.add_argument('--fig2', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 2.')

parser.add_argument('--fig3', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 3.')

parser.add_argument('--fig4', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 4.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')




def combine_svgs(figures, outdir):
    '''Combine svg figs and save.
    '''

    return None


#####MAIN#####
args = parser.parse_args()
fig1 = args.fig1[0]
fig2 = args.fig2[0]
fig3 = args.fig3[0]
fig4 = args.fig4[0]
outdir = args.outdir[0]

Figure("24cm", "6.5cm",
       Panel(
            SVG(fig1).scale(0.5),
            Text("(a)", 35, 20, size=12, weight='bold') #Ttext after if moving/scaling
          ),
       Panel(
          SVG(fig2).scale(0.5),
          Text("(b)", 35, 20, size=12, weight='bold')

        ),Panel(
            SVG(fig3).scale(0.5), #0.318 if 11x11
             Text("(c)", 35, 20, size=12, weight='bold') #Ttext after if moving/scaling
           )).tile(3,1).save("figureS2.svg")
       # Panel(
       #    SVG(fig4).scale(0.35).move(280, 280),
       #    Text("(d)", 25, 20, size=12, weight='bold').move(280,280)
       # )).save("figure2.svg")

pdb.set_trace()
