# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:10:55 2019

@author: NickolasT
"""
import os
from leafSolver import *
import translator as tr
import grapher as gr
import matplotlib.pyplot as plt
import numpy as np
import time
import pydicom
import mapper as mp

globalTime = time.time()
# This is the main file, excute this python script
# Include any necessary function calls here for testing purposes

# =====================PARAMETERS=====================
nLvls = 6                           # Number of layers to use
xStep = 0.25                        # x descretization
searchStep = 1E-2                   # Step size for search values
weightSearchStep = 1E-2                # Threshold of weight convergence
maxWeight = 3                       # Maximum allowed weight
# ==============CONTROLS======================
DO_LEAF_POSITION_SOLVE = True       # True=Run solver.py methods to optimize leaf positions
                                    # (this only needs to be ran once)
DO_LEAF_EXPORT = False              # Use export functions from translator.py
                                    # to export leaf and weight info
rowEx = 11                           # Row index to display 1D fluence and leaf positions
tableWidth = 20                     # Char width of table entries

mlcplan = pydicom.dcmread(r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9FldHNMLC.dcm")

if DO_LEAF_POSITION_SOLVE:
    tempTime = time.time()
    idealFluencesArray = []
    solverAvgArray = []
    #right now seems that it doesnt like the fully blocked locations to have anything but 0 as the value
    
    for f,beam in enumerate(mlcplan.BeamSequence):
        print("Field: " + str(f+1) + " of " + str(len(mlcplan.BeamSequence)))
        fieldTime = time.time()
        shortname = str(beam.BeamNumber)
        out = tr.importDicomFluence(beam)
        xMin=out[5]/10
        xMax=xMin + (out[1]-1)*out[3]/10
        params = (nLvls, xMin, xMax, xStep,searchStep, maxWeight, weightSearchStep)
        # Read in the field data

        optimalFluence = out[-1]
        
        # Simplify the field data by a factor of 4
        # Results in 1cm leaves right now if the fluences from eclipse come over as 0.25cm
        #normFluence = tr.normalizedKernel(optimalFluence,4)
        #use factor of 2 for 0.5cm leaves
        normFluence = tr.normalizedKernel(optimalFluence,2)
        nRows = len(normFluence)
        
        # Instantiate a solver object
        mySolver3 = leafSolver(params)
        solve3, ideal, x, time3 = mySolver3.solveExtend2PeakWeightLimit(normFluence)
        idealFluencesArray.append(ideal)       
        solverAvgArray.append(mySolver3)

    t = (time.time() - tempTime)
    print("*******LEAF SOLVE COMPLETED***********  T: "+ str(t)[0:4] +" s")
    
    if DO_LEAF_EXPORT:    
        for w,mySolve in enumerate(solverAvgArray[0]):
            filenameLeafOut = str(w)+"_"+str(nLvls)+"lvls[FINAL].leaves"       # File to export leaf positions
            mySolve.weight = solverAvgArray[1][w]
            tr.outputLeaves(filenameLeafOut, mySolve,fileDelim)

#now write the data back to the dicom file
tr.convertToComp(mlcplan, solverAvgArray, "HNPreWeightCompPlan.dcm")

print("*******EXPORT COMPLETED***********")

# plt.imshow(idealFluencesArray[0])
# mlmlcfluences=mp.getSolverArrayFluences(solverAvgArray)
# plt.imshow(mlmlcfluences[0])