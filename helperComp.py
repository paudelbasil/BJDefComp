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
import time

from scipy.optimize import minimize
from ansys.mapdl.core import launch_mapdl
from ansys.mapdl import reader as Reader

dbg = False


# Calls APDL and runs the simulation, return result object
def RunSimulation(mapdl, wdir, mainFile):
    '''
    

    Parameters
    ----------
    mapdl : ansys.mapdl
        ansys.mapdl object.
    wdir : string
        working directory for mapdl.
    mainFile : TYPE
        main input file for mapdl.

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    '''
    t=tic()
    ShowMsg(' Simulation started on : ' + getTime())
    jname = 'sJob_%5d' % t
    
    mapdl.clear()
    #mapdl.jobname = jname
    monitor=mapdl.input(FullPath(mainFile,wdir), verbose= False, progress_bar=True)
    mapdl.finish
    out = mapdl.result
    ShowMsg(' Simulation completed on : ' + getTime())
    tok(t)
    
    return out

def CreateInput(wdir, simZeroFile, compFactor, outNodeFile):
    allResult = GetNodalData(FullPath(simZeroFile, wdir))
    
    nodeCnt = len(allResult)
    outputFile = FullPath(outNodeFile, wdir)
    
    f=open(outputFile,"w")  
    
    print("%s %68s" % ('!',"Nodal Information Dump"),file=f)
    print("%s %58s%10i" % ('!',"Total Nodes =", nodeCnt),file=f)
    print("%s%9s%20s%20s%20s" % ('!',"Node","X","Y","Z"),file=f)
    print("nblock,3,,%i" % (nodeCnt),file=f)
    print("(1i9,3e20.9e3)",file=f)

    # Separate node locations and nodal displacements
    nodes = allResult[:,1:4]
    nodeDisp = allResult[:,4:7]    
    for j in range(1,nodeCnt+1):
        print("%9i%20.9E%20.9E%20.9E" % 
              (j, nodes[j-1,0]+compFactor*nodeDisp[j-1,0],
               nodes[j-1,1]+compFactor*nodeDisp[j-1,1],
               nodes[j-1,2]+compFactor*nodeDisp[j-1,2]),file=f)
        
    print("-1",file=f)
    print("! ====================================================================",file=f)
    f.close()
    
### Use the Optimization routine to find the right scale factor instead
### Use scipy.optimize to minimize
# First create a objective function to simplify the iterator
def Objective(inputX, mapdl, wdir, mainFile, simZeroResult, target):
    ShowMsg(".",2)
    ShowMsg("Iteration input Scale factor: %f" % (inputX))
    nodalFileName = 'nodalData_IT.inp'
    nodalFileIt = 'COMP_RSLT_IT_%f.csv' % np.round(inputX,4) 
    
    # Create compensated input file for the prescribed scale
    # SaveNodalInfo(dispZero, nodalFileName, wdir, -inputX, False)
    CreateInput(wdir, simZeroResult, -inputX, nodalFileName)
    
    #print(' -- IN -- ' )
    # inPoints = GetNodalData(FullPath(nodalFileName, wdir))
    # PlotPoints(inPoints, None, 180, -90)

    if(os.path.exists(FullPath(nodalFileIt,wdir))==False):
        # Step: 1 : Run the simulation
        
        # Run the simulation for prescribed scaleFactor
        result = RunSimulation(mapdl, wdir, mainFile)      # Same Result Object as the Reader class
        # Save simulation result
        SaveNodalInfo(result, nodalFileIt, wdir, 1, False) # Save Iteration result
    
        #current = GetNodalData_Result(result, 1.0)
    else:
        # Simulation results already exists, so retrieve from storage
        ShowMsg(' Using data from previously saved file.',3)
    
    # Get current simulation results from file
    current = GetNodalData(FullPath(nodalFileIt,wdir))
        
    # Plot initial and final
    
    # print(' -- OUT -- ' )
    # plt.subplot(1,2,2)
    # outPoints = GetNodalData(FullPath(nodalFileIt, wdir))
    # PlotPoints(outPoints,None, 180, -90)
    
    # Step : 2 : Calculate Error
    nodalError, residualError = CompareNodes(target, current)
        
    return residualError
        
### EXPORT GEOMETRY FROM LAST SOLUTION
def Iterate(mapdl, target, errTol = 1e-5, maxItn=10, wdir='', mainFile =''):
    
    # Do initial simulation
    simZeroResult = 'RESULT_IT0.csv'
    if(os.path.exists(FullPath(simZeroResult,wdir))==False):
        
        rslt = RunSimulation(mapdl, wdir, 'mainInputBox0.dat')
        SaveNodalInfo(rslt, simZeroResult, wdir, readFromFile=False, doCombined=False)
    else:
        ShowMsg(' Using a previously stored BASE SOLUTION.',3)
        
    # Call the scipy minimizer
    initialScale = 1
    
    # DOCS: https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#unconstrained-minimization-of-multivariate-scalar-functions-minimize
    # Use NON-Gradient based iterative algorithm
    # Options: xatol (absolute tolerance in x), method='nelder-mead'
                      
    result = minimize(Objective, initialScale, 
                      (mapdl, wdir, mainFile, simZeroResult, target),
                      tol = errTol, method='nelder-mead',
                      options={'maxiter': maxItn,'disp':True})
    
    # Summarize the result
    ShowMsg('------------------------------------------------------',-1)
    ShowMsg('Status : %s' % result['message'],2)
    ShowMsg('Total Iterations : %d' % result['nfev'],2)
    ShowMsg('Optimum Value    : %d' % result['x'],2)
    ShowMsg('------------------------------------------------------',-1)
    
    # evaluate solution
    # solution = result['x']
    # evaluation = objective(solution)
    #print('Solution: f(%s)' % (solution))
    ShowMsg('Iteration completed. ---')
    
def CompareNodes(target, current):
    '''
    

    Parameters
    ----------
    target : TYPE
        Item 1 to compare.
    current : TYPE
        Item 2 to compare.

    Returns
    -------
    list
        [nodalMax, globalMax].

    '''
    # Both data structure are [N, X, Y, Z]
    nodes = target[:]
    nodeCount  = len(nodes)
    if(dbg):
        ShowMsg('Rows count : %d, Cols count : %d' % (nodeCount,len(target[0])))
    
    nodalDist = nodeCount * [0]
    
    nodalMax, globalMax = 0.0,0.0 
    
    for nodeid in range(nodeCount):
        dx = target[nodeid][1] - current[nodeid][1]
        dy = target[nodeid][2] - current[nodeid][2]
        dz = target[nodeid][3] - current[nodeid][3]
        dist = (dx**2+dy**2+dz**2)**0.5
        nodalDist[nodeid] = dist
        nodalMax = max(nodalMax, dist)
        globalMax = globalMax + dist**2
    
    # Normalize globalMax
    globalMax = globalMax**(0.5)
    globalMax = globalMax
    
    ShowMsg("Max Nodal Diff: %10f Residual Error: %10f" % (nodalMax, globalMax),3)
    
    if(dbg):
        print([np.shape(current), np.shape(target)])
    
    PlotPoints(current, nodalDist,180,-90)
    return [nodalMax, globalMax]

def SaveDefHistory(rstfile, outFile, wdir='', factor=1, readFromFile=True):
    if(readFromFile):
        if(os.path.exists(rstfile)==False):
            ShowMsg('File does NOT exist.')
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
        ShowMsg('Step %i of %i written.' % (ld, rsetMax-1), 1)
    
    # Export the Final solution
    if((rsetMax-ld)>1):
        SaveNodalInfo(rstfile, outFile, wdir,factor,readFromFile, rsetMax-1, True)
        ShowMsg('Step %i of %i written.' % (rsetMax-1, rsetMax-1), 1)
    
    return None
     
def SaveNodalInfo(rstfile, outFile, wdir='', factor=1, readFromFile=True, step=-1, appendTimeStamp=False, doCombined=True):
    
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
    
    # doCombined = True
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
    fullname= FullPath(outFile, wdir) 
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
    '''
    

    Parameters
    ----------
    result : mapdl.result
        Result object from mapdl.
    factor : float, optional
        Scaling for displacement. The default is 1.0.

    Returns
    -------
    nodalData : array(:,4)
        Nodal data from the mapdl.result with scaled displacments.

    '''
    rsetMax=result.nsets-1
    
    # Plot nodal solution
    # result.plot_nodal_solution(rsetMax,background='w',show_edges=True,show_displacement=True)
    
    # Get nodes original position
    nodes = result.mesh.nodes
    nodeCnt = len(nodes)
    
    compFactor = factor       # To get deformed shape, use +1.
    rset = rsetMax
    # tStep = result.solution_info(rset)['timfrq']
    
    # Get nodal solutions, along with nodal numbers
    nodeNum, nodeDisp = result.nodal_displacement(rset)   # last set
    nodalData = nodeCnt*[[0,0,0,0]]
    for j in range(nodeCnt):
        nX = nodes[j-1,0] + compFactor*nodeDisp[j-1,0]
        nY = nodes[j-1,1] + compFactor*nodeDisp[j-1,1]
        nZ = nodes[j-1,2] + compFactor*nodeDisp[j-1,2]
        nodalData[j]= [j+1, nX, nY, nZ]
        
    # Return 2D array instead of list of list
    out = np.array(nodalData)
    return out

def PlotNodalSolution(rstfile, step=-1, readFromFile=True):
    '''
    

    Parameters
    ----------
    rstfile : TYPE
        Either filepath of 'rst' or the result object variable.
    step : TYPE, optional
        DESCRIPTION. The default is -1.
    readFromFile : bool, optional
        If set to False, rstfile is 'mapdl.result' object. The default is True.

    Returns
    -------
    result : PyVista.pyplot
        Plot of nodal deformation.

    '''
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
    result.plot_nodal_solution(rsetMax,background='w',
                               show_edges=True,
                               show_displacement=True, 
                               notebook=False, 
                               window_size=[512,384])
    return result

def PlotPoints2(nodes):
    
    ptCloud = nodes[:,1:4]  # Disregard the first column (nodeIDs)
    pc = pv.PolyData(ptCloud)
    pc.plot(background='w',notebook=False,render_points_as_spheres=True,
            eye_dome_lighting=True, window_size=[512,384], parallel_projection=True)
    return pc

def GetNodalData(file):
    '''
    

    Parameters
    ----------
    file : str
        Full file path of the csv file.

    Returns
    -------
    dfarray : array(array())
        Table of Nodal location [[node, x, y, z],].

    '''
    # import pandas as pd
    df = pd.read_csv(file, engine="python", skiprows=5, skipfooter=2, 
                     comment='!', sep='\s+', skipinitialspace=True, 
                     header=None)
    dfarray = df.values
    return dfarray

def FullPath(file, wdir):
    '''
    

    Parameters
    ----------
    file : str
        Full file path (excluding dir).
    wdir : TYPE
        Directory with/without closing backslash.

    Returns
    -------
    TYPE
        Full file path (with dir).

    '''
    if(wdir[-1] != "\\"):
        wdir = wdir + "\\"
    return wdir + file

def SplitFileName(fname):
    '''
    

    Parameters
    ----------
    fname : Str
        filename (without directory).

    Returns
    -------
    list
        [name, ext].

    '''
    fileExt = fname.split(".")[-1]
    fileName = fname.split('.'+fileExt)[0]
    return [fileName, fileExt]

def PlotPoints(point3d, cdata=None, v1=0,v2=0, ax = None):
    """    

    Parameters
    ----------
    point3d : TYPE
        [[id, x, y, z],...].
    v1 : TYPE, optional
        Viewer - X rotation . The default is 0.
    v2 : TYPE, optional
        Viewer - Y rotation. The default is 0.

    Returns
    -------
    Static Matplotlib Plot .

    """
    
    if(ax==None):
        ax = plt.axes(projection='3d')
    
    # 0 : Node Number
    xdata = point3d[:,1]
    ydata = point3d[:,2]
    zdata = point3d[:,3]
    
    if(cdata==None):
        cdata = xdata+ydata+zdata
    
    pp=ax.scatter3D(xdata, ydata, zdata, c=cdata, cmap='jet', marker=',')
    ax.view_init(v1, v2)
    ax.grid(False)
    ax.set_proj_type('ortho')
    plt.axis('off')
    plt.colorbar(pp)
    plt.show()
    return ax

def tic():
    return time.time()

def tok(t):
    el = tic() - t
    ShowMsg(' Elapsed time : %8.3f s' % el)
    return

def getTime():
    ltime = time.localtime(time.time())
    return time.asctime(ltime)
    
def ShowMsg(msg='', level=0):
    if(level == 0):
        print('>>> ' + str(msg))
    elif(level == 1):
        print('--- ' + str(msg))
    elif(level == 2):
        print('\t ' + str(msg))
    elif(level == 3):
        print('*** ' + str(msg))
    else:
        print(' -----------------------------------------------------')
    
    return None

