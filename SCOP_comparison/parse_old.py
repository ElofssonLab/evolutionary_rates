#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import os
import argparse
from scipy import ndimage

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses the 2009 alignments and returns
                                                the SCOP identifiers as well as alignments for tree-puzzle.''')

parser.add_argument('--infile', nargs=1, type= str, default=sys.stdin, help = 'path to infile.')



#####MAIN#####
args = parser.parse_args()
infile = args.infile[0]

parse(infile)
