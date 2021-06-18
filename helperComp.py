# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 21:12:39 2021

@author: bjp89
"""
# Imports
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from ansys.mapdl.core import launch_mapdl
from ansys.mapdl import reader as Reader


### EXPORT GEOMETRY FROM LAST SOLUTION
def Iterate(mapdl, target, errTol = 1e-5, maxItn=10, wdir='', mainFile =''):
    itn = 0 
    residualError, nodalError    = 1e6,1e6
    
    ## ITERATION LOOP
    while(residualError >errTol and nodalError > errTol and itn<maxItn):
        print(">>> Iteration loop : %i" % itn)
        nodalFileIt = 'COMP_RSLT_IT_%i.csv' % itn 
        # Run the simulation
        mapdl.clear()
        mapdl.input(wdir+mainFile, verbose=False)
        mapdl.finish()
        result = mapdl.result       # Same Result Object as the Reader class
        # print(result)  # Display result for futher info
        
        itn = itn + 1
        
        # current = GetNodalData(nodalFileIt)
        current = GetNodalData_Result(result, 1.0)
        # Check Error
        nodalError, residualError = CompareNodes(target, current)
        
        # Compensate geom
        SaveNodalInfo(result, nodalFileIt, wdir, 1, False) # Save History
        SaveNodalInfo(result, 'nodalData.inp', wdir, -0.95, False)
    
    # Check if converged
    if(itn>=maxItn):
        # Not Converged
        print('--- Compensation algorithm did NOT converge')
    else:
        # Converged
        print('--- Compensation algorithm converged.')
        SaveNodalInfo(result, 'nodalData.inp', wdir, 0, False)
    
    print('>>> Iteration completed. ---')

    
def CompareNodes(target, current):
    # Both data structure are [N, X, Y, Z]
    nodes = target[0]
    nodeCount  = len(nodes)
    nodalDist = nodeCount * [0]
    
    nodalMax, globalMax = 0,0 
    
    for nodeid in range(nodeCount):
        dx = target[1][nodeid] - current[1][nodeid]
        dy = target[2][nodeid] - current[2][nodeid]
        dz = target[3][nodeid] - current[3][nodeid]
        dist = (dx**2+dy**2+dz**2)**0.5
        nodalDist[nodeid-1] = dist
        nodalMax = max(nodalMax, dist)
        globalMax = globalMax + dist**2
    
    # Normalize globalMax
    globalMax = globalMax**(0.5)
    globalMax = globalMax/nodeCount
    
    print("Max Nodal Diff: %10f Residual Error: %10f" % (nodalMax, globalMax))
    return [nodalMax, globalMax]

def SaveDefHistory(rstfile, outFile, wdir='', factor=1, readFromFile=True):
    if(readFromFile):
        if(os.path.exists(rstfile)==False):
            print('File does NOT exist.')
            return None
        # Create result object by loading the result file
        # Working with result after simulation
        result = Reader.read_binary(rstfile)  # Using Reader
    else:
        result = rstfile

    rsetMax = result.nsets

    loads = range(0,rsetMax-1,10)
    
    for ld in loads:
        SaveNodalInfo(rstfile, outFile, wdir,factor,readFromFile, ld, True)
        print('--- Step %i of %i written.' % (ld, rsetMax-1))
    
    # Export the Final solution
    if((rsetMax-ld)>1):
        SaveNodalInfo(rstfile, outFile, wdir,factor,readFromFile, rsetMax-1, True)
        print('--- Step %i of %i written.' % (rsetMax-1, rsetMax-1))
    
    return None
     

def SaveNodalInfo(rstfile, outFile, wdir='', factor=1, readFromFile=True, step=-1, appendTimeStamp=False):
    
    if(readFromFile):
        # Create result object by loading the result file
        # Working with result after simulation
        result = Reader.read_binary(rstfile)  # Using Reader
    else:
        result = rstfile

    if(step<0):    
        rsetMax=result.nsets-1
    else:
        rsetMax = step
    
    # Get nodal displacements
    nnum,ndisp = result.nodal_displacement(rsetMax)
    
    # Plot nodal solution
    # result.plot_nodal_solution(rsetMax,background='w',show_edges=True,show_displacement=True)
    
    # Get nodes original position
    nodes = result.mesh.nodes
    
    nodeCnt = len(nodes)
    
    doCombined = True
    compFactor = factor       # To get deformed shape, use +1.
    rset = rsetMax
    tStep = result.solution_info(rset)['timfrq']
    
    tStepStr = "%08.1f" % (tStep)
    
    # Get nodal solutions, along with nodal numbers
    nodeNum, nodeDisp = result.nodal_displacement(rset)   # first set
    nodeNum, nodeTemp = result.nodal_temperature(rset)
    isothermTemp = nodeTemp[1]                          # Nodal temp
    
    if (appendTimeStamp):
        filePref, fileExt = SplitFileName(outFile)
        outFile = filePref + '_T_' + tStepStr + '.' + fileExt

    # Using info https://www.w3schools.com/python/python_file_write.asp 
    fullname= wdir + outFile
    f=open(fullname,"w")  
    
    print("%s %68s" % ('!',"Nodal Information Dump"),file=f)
    print("%s %58s%10i" % ('!',"Total Nodes =", nodeCnt),file=f)
    if(doCombined):
        print("%s%9s%20s%20s%20s" % ('!',"Node","X","Y","Z"),file=f)
        print("nblock,3,,%i" % (nodeCnt),file=f)
        print("(1i9,3e20.9e3)",file=f)
    
    else:
        print("%s%9s%20s%20s%20s%20s%20s%20s" % ('!',"Node","X","Y","Z","UX","UY","UZ"),file=f)
        print("nblock,3,,%i" % (nodeCnt),file=f)
        print("(1i9,6e20.9e3)",file=f)
        
    for j in nodeNum:
        if(doCombined):
            print("%9i%20.9E%20.9E%20.9E" % (j,nodes[j-1,0]+compFactor*nodeDisp[j-1,0],nodes[j-1,1]+compFactor*nodeDisp[j-1,1],nodes[j-1,2]+compFactor*nodeDisp[j-1,2]),file=f)
        else:
            print("%9i%20.9E%20.9E%20.9E%20.9E%20.9E%20.9E" % (j,nodes[j-1,0],nodes[j-1,1],nodes[j-1,2],nodeDisp[j-1,0],nodeDisp[j-1,1],nodeDisp[j-1,2]),file=f)
        
    print("-1",file=f)
    print("! ====================================================================",file=f)
    f.close()
    return True

def GetNodalData_Result(result, factor=1.0):
    rsetMax=result.nsets-1
    
    # Get nodal displacements
    nnum,ndisp = result.nodal_displacement(rsetMax)
    
    # Plot nodal solution
    # result.plot_nodal_solution(rsetMax,background='w',show_edges=True,show_displacement=True)
    
    # Get nodes original position
    nodes = result.mesh.nodes
    nodeCnt = len(nodes)
    
    compFactor = factor       # To get deformed shape, use +1.
    rset = rsetMax
    # tStep = result.solution_info(rset)['timfrq']
    
    # Get nodal solutions, along with nodal numbers
    nodeNum, nodeDisp = result.nodal_displacement(rset)   # first set
    nodalData = nodeCnt*[[0,0,0,0]]
    for j in range(nodeCnt):
        nodalData[j]= [j+1, nodes[j-1,0]+compFactor*nodeDisp[j-1,0], nodes[j-1,1]+compFactor*nodeDisp[j-1,1], nodes[j-1,2]+compFactor*nodeDisp[j-1,2] ]
        
    return nodalData

def PlotNodalSolution(rstfile, step=-1, readFromFile=True):
    if(readFromFile):
        # Create result object by loading the result file
        # Working with result after simulation
        result = Reader.read_binary(rstfile)  # Using Reader
    else:
        result = rstfile

    if(step<0):
        rsetMax=result.nsets-1
    else:
        rsetMax = step
        
    # Plot nodal solution
    result.plot_nodal_solution(rsetMax,background='w',show_edges=True,show_displacement=True, notebook=False, window_size=[512,384])
    return result

def GetNodalData(file):
    # import pandas as pd
    df = pd.read_csv(file, engine="python", skiprows=5, skipfooter=1, comment='!', sep='\s+', skipinitialspace=True, header=None)
    dfarray = df.values
    return dfarray

def FullPath(file, wdir):
    return wdir + file

def SplitFileName(fname):
    fileExt = fname.split(".")[-1]
    fileName = fname.split('.'+fileExt)[0]
    return [fileName, fileExt]


def PlotPoints(point3d, v1=0,v2=0):
    ax = plt.axes(projection='3d')
    
    # 0 : Node Number
    xdata = point3d[:,1]
    ydata = point3d[:,2]
    zdata = point3d[:,3]
    cdata = xdata+ydata+zdata
    ax.scatter3D(xdata, ydata, zdata, c=cdata, cmap='inferno')
    ax.view_init(v1, v2)
    ax.grid(False)
    
    