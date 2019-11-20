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

parser.add_argument('--fig5', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 5.')

parser.add_argument('--fig6', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 6.')

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
fig5 = args.fig5[0]
fig6 = args.fig6[0]
outdir = args.outdir[0]

Figure("21cm", "16.5cm",
       Panel(
            SVG(fig1).scale(0.35).move(10,0),
            Text("(a)", 25, 20, size=12, weight='bold') #Ttext after if moving/scaling
          ),
       Panel(
          SVG(fig2).scale(0.35).move(10,0),
          Text("(b)", 25, 20, size=12, weight='bold')
          ),
        Panel(
             SVG(fig3).scale(0.35).move(10,0), #0.318 if 11x11
              Text("(c)", 25, 20, size=12, weight='bold') #Ttext after if moving/scaling
            ),
        Panel(
           SVG(fig4).scale(0.35).move(10,0),
           Text("(d)", 25, 20, size=12, weight='bold')
        ),
        Panel(
             SVG(fig5).scale(0.35).move(10,0), #0.318 if 11x11
              Text("(e)", 25, 20, size=12, weight='bold') #Ttext after if moving/scaling
            ),
        Panel(
           SVG(fig6).scale(0.35).move(10,0),
           Text("(f)", 25, 20, size=12, weight='bold')
           )).tile(3,2).save("figureS2.svg")
pdb.set_trace()
