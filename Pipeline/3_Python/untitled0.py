# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:22:56 2021

@author: rmoge
"""

from __future__ import division
import re, os,sys, math, operator,random, copy,collections,time; import numpy as np; import pandas as pd; import csv
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as smt

from pathlib import Path
homedir = str(Path.home())
workingdir = homedir + '\\Documents\\GitHub\\Bifidobacterium'
os.chdir(workingdir)
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

#%%
import goatools as goat

from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo

fin_dag = download_go_basic_obo("go-basic.obo")

godag = GODag(fin_dag)

goid = 'GO:0050807'


def prt_flds(gosubdag):
    """Print the available printing fields"""
    print('Print fields:')
    for fld in sorted(gosubdag.prt_attr['flds']):
        print('    {F}'.format(F=fld))
        
        

from goatools.gosubdag.gosubdag import GoSubDag

# Create a subset of the GO DAG which contains:
#   * The selected GO term and
#   * All the GO terms above it
gosubdag = GoSubDag(goid, godag, relationships=True, prt=False)

# Get additional information for chosen GO
ntgo = gosubdag.go2nt[goid]

# Choose fields and custom printing format
# prt_flds(gosubdag)  # Uncomment to see the available print fields
prtfmt = '{NS} {GO} D{depth:02} {GO_name}'

# Print detailed information for GO
print(prtfmt.format(**ntgo._asdict()))        