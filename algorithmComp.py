# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 17:46:29 2021

@author: bjp89
"""
# Imports
import os
import shutil

from helperComp import * # Helper function for Iterations

from ansys.mapdl.core import launch_mapdl
from ansys.mapdl import reader as Reader

# This will server as the template to run Compensation algorithm for use 
#       with BJ Sintering Simulation

# Define all parameters

mainFile = 'mainInputBox.dat'
nodalFile = 'nodalData.inp'
wdir = 'D:\\Basil\\Deformation_Comp\\Hao\\CD4\\'
nodalFileIt = ''
rstFile = ''

# Backup Original Nodal Information
nodalNameOnly = os.path.splitext(nodalFile)[0]
nodalExtOnly = os.path.splitext(nodalFile)[1]
backupFile = nodalNameOnly + "_ORIG" + nodalExtOnly
backupNodalInfo = shutil.copyfile(wdir +  nodalFile,  wdir +  backupFile)

# STORE ORIGINAL GEOMETRY CONFIGURATION AS TARGET
targetNodes = GetNodalData( wdir + backupFile)

# Load the solver's server
mapdl = launch_mapdl(run_location=wdir,jobname='file', nproc=2)
# rstFile = FullPath('file.rst', wdir)

# Run Iteration,
#  Iterate(mapdl, target, errTol = 1e-5, maxItn=10, wdir='', mainFile ='')
Iterate(mapdl, targetNodes,1e-5,5,wdir, mainFile)
    
mapdl.exit()    # Close the APDL
    
## EXPORT DEFORMED GEOMETRY

## COMPARE SOLUTION

## CONVERGED SOLUTION 
### EXPORT PRE-DEFORMED GEOMETRY

## UN-CONVERGED SOLUTION
