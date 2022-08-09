# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:10:55 2019

@author: NickolasT and Andrew Morrow
"""
from leafSolver import *
import translator as tr
import numpy as np
import time
import pydicom

# This is the main file, excute this python script
# Include any necessary function calls here for testing purposes

# =====================PARAMETERS=====================
nLvls = 2                           # Number of layers to use
xStep = 0.25                        # x discretization.  use 0.25cm when doing IMRT mode cuz that's what the fluence is pieced up as
searchStep = 1E-2                   # Step size for search values
weightSearchStep = 1E-2                # Threshold of weight convergence
maxWeight = 3                       # Maximum allowed weight
plantype = "ra"

if plantype == "ra":
    dlg = 0.1103 #the killeen 6x DLG is 0.1103. 
    dlg=dlg/2 #per leaf
    raplan = pydicom.dcmread(r"K:\Physics Division\Personal folders\ANM\MLMLC\data\fortesting\idealplan.dcm")
    tempTime = time.time()
    idealFluencesArray = []
    solverAvgArray = []
    #for rapid arc plans feeding the algorithm, we have to sum together the controlpoints across the fields in the plan when they are at the same gantry angle.
    #then we need to fit to that.
    #we'll have to make a new RT dicom plan as well at the end.
    #we are not allowing for jaw tracking right now, so xmin and xmax are in the first control point
    xMin=raplan.BeamSequence[0].ControlPointSequence[0].BeamLimitingDevicePositionSequence[0].LeafJawPositions[0]
    xMax=raplan.BeamSequence[0].ControlPointSequence[0].BeamLimitingDevicePositionSequence[0].LeafJawPositions[1]
    yMin=raplan.BeamSequence[0].ControlPointSequence[0].BeamLimitingDevicePositionSequence[1].LeafJawPositions[0]
    yMax=raplan.BeamSequence[0].ControlPointSequence[0].BeamLimitingDevicePositionSequence[1].LeafJawPositions[1]
    #ypos=raplan.BeamSequence[0].BeamLimitingDeviceSequence[2].LeafPositionBoundaries #positions of the leaves in the Y direction, or just hard code it in.
    #hard codingg it in.  dont have fields bigger than 15cm (x) x 20cm(y) for now.  y goes between MLC indices 10 and 50 for half cm leaves.
    mus=[]
    nx = np.int16(np.round(abs(np.round((xMax+1-xMin)/xStep))))+1
    ny = np.int16(np.round(abs(np.round((yMax+1-yMin)/5))))
    
    for b in range(len(raplan.FractionGroupSequence[0].ReferencedBeamSequence)):
        mus.append(raplan.FractionGroupSequence[0].ReferencedBeamSequence[b].BeamMeterset)
    for cp, controlpoint in enumerate(raplan.BeamSequence[0].ControlPointSequence):
        if cp == 0:
            continue#skip first cp - it just has jaw positions in it at least with no jaw tracking
        print("CP: " + str(cp) + " of " + str(len(raplan.BeamSequence[0].ControlPointSequence)))
        farray=np.zeros((ny, nx),dtype=np.float64)
        for f,beam in enumerate(raplan.BeamSequence):

            leaves = raplan.BeamSequence[f].ControlPointSequence[cp].BeamLimitingDevicePositionSequence[0].LeafJawPositions
            #weight will be the mu/cp
            cpmu = np.float64(mus[f]*(raplan.BeamSequence[0].ControlPointSequence[cp].CumulativeMetersetWeight-raplan.BeamSequence[0].ControlPointSequence[cp-1].CumulativeMetersetWeight))
            farray = np.add(farray,tr.convertCPtoFluence(leaves,cpmu,xStep,xMin, xMax, yMin, yMax, dlg))
            mus.append(cpmu)
        
        params = (nLvls, xMin, xMax, xStep,searchStep, maxWeight, weightSearchStep)
        normFluence = tr.normalizedKernel(farray,1)
        
        # Instantiate a solver object
        mySolver3 = leafSolver(params)
        solve3, ideal, x, time3 = mySolver3.solveExtend2PeakWeightLimit(normFluence)     
        solverAvgArray.append(mySolver3)
    t = (time.time() - tempTime)
    print("*******LEAF SOLVE COMPLETED***********  T: "+ str(t)[0:4] +" s")
    #now write the data back to the dicom file
    #tr.convertToComp(mlcplan, solverAvgArray, "HNPreWeightCompPlan.dcm")  #need to make a custom function here

if plantype == "imrt":
    mlcplan = pydicom.dcmread(r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9fldhnkideal\idealplan.dcm")
    tempTime = time.time()
    idealFluencesArray = []
    solverAvgArray = []
    #right now seems that it doesnt like the fully blocked locations to have anything but 0 as the value
    for f,beam in enumerate(mlcplan.BeamSequence):
        print("Field: " + str(f+1) + " of " + str(len(mlcplan.BeamSequence)))
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

    #now write the data back to the dicom file
    tr.convertToComp(mlcplan, solverAvgArray, "HNPreWeightCompPlan.dcm")

# plt.imshow(idealFluencesArray[0])
# mlmlcfluences=mp.getSolverArrayFluences(solverAvgArray)
# plt.imshow(mlmlcfluences[0])