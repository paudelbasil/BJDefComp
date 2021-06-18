# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 11:19:00 2021

@author: Basil J. Paudel
"""
from helperComp import *
import os

# --------------------------------------------------------------------------
# Snippet for EXTRACTING results and store compensated geometry for various 
#   compensation factors
# --------------------------------------------------------------------------

mainDir = 'E:\Basil Paudel\Scratch\CD4_V2\\'
# First Data Extraction
wdir = mainDir + 'CD4_V2_IT0_Results\\' 
if(os.path.exists(wdir)==False): os.mkdir(wdir)

rstFile = mainDir + 'CD4_V2_IT0\CD4V2_IT0.rst'
outFile = 'CD4V2_IT0.csv'

SaveDefHistory(rstFile, outFile, wdir)

# Iteration Data Extraction
for i in range(1,4):
    wdir = mainDir + 'CD4_V2_IT1.%i_Results\\' % i    
    if(os.path.exists(wdir)==False): os.mkdir(wdir)
    
    rstFile = mainDir + 'CD4_V2_IT1.%i\CD4V2_IT1.rst' % i
    outFile = 'CD4V2_IT1.%i.csv' % i
    
    SaveDefHistory(rstFile, outFile, wdir)

# --------------------------------------------------------------------------
# Snippet for Reading Result file and store compensated geometry for various 
#   compensation factors
# --------------------------------------------------------------------------

# Results from Iteration 0
    
rstfile = 'E:\Basil Paudel\Scratch\CD4_V2\CD4_V2_IT0\CD4V2_IT0.rst'
outFile = 'nodes_it0.inp'

# Create Nodes with 0.5 scale of actual distortion in IT0
wdir = 'E:\Basil Paudel\Scratch\CD4_V2\CD4_V2_IT1.1\\'
SaveNodalInfo(rstfile, outFile, wdir, factor=-0.5)

# Create Nodes with 1.0 scale of actual distortion in IT0
wdir = 'E:\Basil Paudel\Scratch\CD4_V2\CD4_V2_IT1.2\\'
SaveNodalInfo(rstfile, outFile, wdir, factor=-1)

# Create Nodes with 1.5 scale of actual distortion in IT0
wdir = 'E:\Basil Paudel\Scratch\CD4_V2\CD4_V2_IT1.3\\'
SaveNodalInfo(rstfile, outFile, wdir, factor=-1.5)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
