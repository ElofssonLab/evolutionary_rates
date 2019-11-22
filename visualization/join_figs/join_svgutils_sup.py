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

parser.add_argument('--fig7', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 7.')

parser.add_argument('--fig8', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 8.')

parser.add_argument('--fig9', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 9.')

parser.add_argument('--fig10', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 10.')

parser.add_argument('--fig11', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 11.')

parser.add_argument('--fig12', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 12.')

parser.add_argument('--fig13', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 13.')

parser.add_argument('--fig14', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 14.')

parser.add_argument('--fig15', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 15.')

parser.add_argument('--fig16', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 16.')

parser.add_argument('--fig17', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 17.')

parser.add_argument('--fig18', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 18.')

parser.add_argument('--fig19', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 19.')

parser.add_argument('--fig20', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 20.')

parser.add_argument('--fig21', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 21.')

parser.add_argument('--fig22', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 22.')

parser.add_argument('--fig23', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 23.')

parser.add_argument('--fig24', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 24.')

parser.add_argument('--fig25', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 25.')

parser.add_argument('--fig26', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 26.')

parser.add_argument('--fig27', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 27.')

parser.add_argument('--fig28', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 28.')

parser.add_argument('--fig29', nargs=1, type= str,
default=sys.stdin, help = 'path to figure 29.')

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
fig7 = args.fig7[0]
fig8 = args.fig8[0]
fig9 = args.fig9[0]
fig10 = args.fig10[0]
fig11 = args.fig11[0]
fig12 = args.fig12[0]
fig13 = args.fig13[0]
fig14 = args.fig14[0]
fig15 = args.fig15[0]
fig16 = args.fig16[0]
fig17 = args.fig17[0]
fig18 = args.fig18[0]
fig19 = args.fig19[0]
fig20 = args.fig20[0]
fig21 = args.fig21[0]
fig22 = args.fig22[0]
fig23 = args.fig23[0]
fig24 = args.fig24[0]
fig25 = args.fig25[0]
fig26 = args.fig26[0]
fig27 = args.fig27[0]
fig28 = args.fig28[0]
fig29 = args.fig29[0]

outdir = args.outdir[0]

#S1
fig = Figure("21cm", "10cm",
        Panel(
             SVG(fig1).scale(0.5), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig2).scale(0.5),
             Text("(b)", 25, 20, size=14, weight='bold')
            )).tile(2,1).save("figureS1.svg")
#S2
#Write text after if moving/scaling
fig = Figure("21cm", "14cm",
        Panel(
             SVG(fig3).scale(0.35), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig4).scale(0.35),
             Text("(b)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig5).scale(0.35), #0.318 if 11x11
             Text("(c)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig6).scale(0.35), #0.318 if 11x11
             Text("(d)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig7).scale(0.35), #0.318 if 11x11
             Text("(e)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig8).scale(0.35), #0.318 if 11x11
             Text("(f)", 25, 20, size=14, weight='bold')
            )).tile(3,2).save("figureS2.svg")

#S3
#Write text after if moving/scaling
fig = Figure("21cm", "14cm",
        Panel(
             SVG(fig9).scale(0.35), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig10).scale(0.35),
             Text("(b)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig11).scale(0.35), #0.318 if 11x11
             Text("(c)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig12).scale(0.35), #0.318 if 11x11
             Text("(d)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig13).scale(0.35), #0.318 if 11x11
             Text("(e)", 25, 20, size=14, weight='bold')
            )).tile(3,2).save("figureS3.svg")

#S4
#Write text after if moving/scaling
fig = Figure("21cm", "14cm",
        Panel(
             SVG(fig14).scale(0.35), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig15).scale(0.35),
             Text("(b)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig16).scale(0.35), #0.318 if 11x11
             Text("(c)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig17).scale(0.35), #0.318 if 11x11
             Text("(d)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig18).scale(0.35), #0.318 if 11x11
             Text("(e)", 25, 20, size=14, weight='bold')
            )).tile(3,2).save("figureS4.svg")

#S5
fig = Figure("21cm", "10cm",
        Panel(
             SVG(fig19).scale(0.5), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig20).scale(0.5),
             Text("(b)", 25, 20, size=14, weight='bold')
            )).tile(2,1).save("figureS5.svg")

#S6
fig = Figure("21cm", "21cm",
        Panel(
             SVG(fig21).scale(0.5), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig22).scale(0.5),
             Text("(b)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig23).scale(0.5), #0.318 if 11x11
             Text("(c)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig24).scale(0.5), #0.318 if 11x11
             Text("(d)", 25, 20, size=14, weight='bold')
            )).tile(2,2).save("figureS6.svg")

#S8
#Write text after if moving/scaling
fig = Figure("21cm", "10cm",
        Panel(
             SVG(fig25).scale(0.4), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig26).scale(0.4),
             Text("(b)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig27).scale(0.4), #0.318 if 11x11
             Text("(c)", 25, 20, size=14, weight='bold')
             )).tile(3,1).save("figureS8.svg")

#S10
fig = Figure("21cm", "10cm",
        Panel(
             SVG(fig28).scale(0.5), #0.318 if 11x11
             Text("(a)", 25, 20, size=14, weight='bold')
            ),
        Panel(
             SVG(fig29).scale(0.5),
             Text("(b)", 25, 20, size=14, weight='bold')
            )).tile(2,1).save("figureS10.svg")
