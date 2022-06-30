# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 23:29:03 2019

@author: NickolasT
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit

# This file calculates the 3D fluence distributions from an inputted fluence projection
# It also returns summed errors when comparing fluence distributions

@jit(nopython=False)
def importTmrMap(filename):
    # Import Tmr map array describing basic beam characteristics as a 2D array
    fd = open(filename,'r')
    lines = fd.readlines()
    xN = len(lines[0].split())
    yN = len(lines)
    
    outMap = np.zeros((yN,xN))
    
    for i,line in enumerate(lines):
        for j,entry in enumerate(line.split()):
            outMap[i][j] = float(entry)
    return outMap, xN, yN

def createCircleMap(radius, xRange,yRange,xN, yN):
    # Creates a 2D array of a map used for a circle mask, where cells inside
    # the circle are 1 and cells outside the circles are 0
    # This is used to creat a cross-section mask for a cylinder volume
    yStep = (yRange[1]-yRange[0])/yN
    xStep = (xRange[1]-xRange[0])/xN
    X = np.arange(xRange[0],xRange[1],xStep)
    Y = np.arange(yRange[0],yRange[1],yStep)
    
    outMap = np.zeros((yN, xN))
    
    for i in range(yN):
        for j in range(xN):
            if (np.sqrt(Y[i]**2 + X[j]**2) <= radius):
                outMap[i][j] = 1
    return outMap, X, Y
@jit(nopython=True)
def transmit1DFluence(inFluence, TmrMap):
    # Takes a 1D array of fluence values (a certain row) and extends it into
    # a 2D array. The extension duplicates the 1D array until it makes a 2D
    # array that has same length of TmrMap. The TmrMap transmission map
    # is then applied.
    out = np.zeros((len(TmrMap),len(inFluence)))
    for i in range(len(TmrMap)):
        out[i] = inFluence
    outMap = out * TmrMap
    return outMap
    
@jit(nopython=True)
def rotFluence(inFluence,xRange, yRange, steps, angle):
    # Takes a 2D array of fluence penetration and rotates the values about the
    # origin (0,0) by the amount of angle (in radians). xRange and yRange are tuples
    # that define the coordinates used, where xRange[0] and yRange[0] are located
    # at the top left of the 2D array. The rotation begins with 0 degrees pointing
    # up, and rotates clockwise. 
    # Rotation of the output values finds the closest corresponding cell from inFluence
    # If a cell is rotated outside of original dimensions, then ignore it
    # If an output cell has no inFluence value rotated into, then set to 0
    yN = len(inFluence)
    if (yRange[1] > yRange[0]):
        yStep = steps[1]
    else:
        yStep = -steps[1]   
    xN = len(inFluence[0])
    if (xRange[1] > xRange[0]):
        xStep = steps[0]
    else:
        xStep = -steps[0]
    transMat = np.array([[np.cos(angle), -np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    X = np.arange(xRange[0],xRange[1]+xStep,xStep)
    Y = np.arange(yRange[0],yRange[1]+yStep,yStep)
    outFluence = np.zeros((yN,xN))
    for i in range(yN):
        for j in range(xN):
            oldCoord = np.array([[X[j]],[Y[i]]])
            newCoord = transMat @ oldCoord
            if (newCoord[0][0] > max(xRange) or newCoord[0][0] < min(xRange)) or (newCoord[1][0] > max(yRange) or newCoord[1][0] < min(yRange)):
                continue
            xIndex = xN-1
            yIndex = yN-1
            for k in range(xN):
                if (newCoord[0][0] < (X[k] + abs(xStep)/2)) and (newCoord[0][0] > (X[k]-abs(xStep)/2)):
                    
                    xIndex = k
            for k in range(yN):
                if (newCoord[1][0] < (Y[k] + abs(yStep)/2)) and (newCoord[1][0] > (Y[k]-abs(yStep)/2)):
                    yIndex = k
            outFluence[i][j] = inFluence[yIndex][xIndex]
            
    return outFluence,xRange,yRange

@jit(nopython=True)
def mapFluences(fluencesArray,solverParams):
    # Takes array of fluences (with indices [angle][row][x]) and returns the
    # fluence distribution in a volume defined by solverParams
    xMin = solverParams[0]
    xMax = solverParams[1]
    xStep = solverParams[2]
    ySize = solverParams[3]
    xSize = solverParams[4]
    maskMap = solverParams[5]
    angleArray = solverParams[6]
    tmrMap = solverParams[7]
    volume = np.zeros((len(fluencesArray[0]), ySize, xSize))
    for a,angle in enumerate(angleArray):
        fluence = fluencesArray[a]
        for i in range(len(fluencesArray[0])):
            pen = transmit1DFluence(fluence[i], tmrMap)
            pen = rotFluence(pen,(xMin,xMax),(xMax,xMin),(xStep,xStep),angle/180*np.pi)[0]
            pen = maskMap * pen
            volume[i] = volume[i] + pen
    return volume

def getSolverArrayFluences(solverArray):
    # Given by array of solvers, returns array of fluences (with indices [angle][row][x])
    # Normalizes the fluence by getting the fluence from the solver with weight of 1
    #nick's code doesnt allow for varying field size
    #fluencesArray = np.zeros((len(solverArray),len(solverArray[0].getFluence(1)),len(solverArray[0].getFluence(1)[0])))
    fluencesArray=[]
    for i in range(len(solverArray)):
        #fluencesArray[i] = solverArray[i].getFluence(1)
        fluencesArray.append(solverArray[i].getFluence(1))
    return fluencesArray
def getSolverArrayFluencesWithWeights(solverArray, weightArray):
    # Given by array of solvers, returns array of fluences (with indices [angle][row][x])
    # The fluence gotten from the solver uses inputted weight array (with indices [angle])
    fluencesArray = np.zeros((len(solverArray),len(solverArray[0].getFluence(1)),len(solverArray[0].getFluence(1)[0])))
    for i in range(len(solverArray)):
        fluencesArray[i] = solverArray[i].getFluence(weightArray[i])
    return fluencesArray
@jit(nopython=True)
def diffVolume(volume1,volume2):
    # Calculates the difference between corresponding fluence volumes of volume1 and volume2
    # Uses the difference to create RGB values
    # Where white is lowest difference, red is volume2 > volume1, and blue is volume2 < volume1
    # Also returns the total error and over error
    sumError = 0
    sumErrorOver = 0
    out = np.zeros((len(volume1), len(volume1[0]), len(volume1[0][0]), 3), dtype=np.uint8)
    diffMax = np.amax(np.abs(volume2-volume1))
    for i in range(len(volume1)):
        for j in range(len(volume1[0])):
            for k in range(len(volume1[0][0])):
                diff = abs(volume2[i][j][k] - volume1[i][j][k])
                if (volume2[i][j][k] > volume1[i][j][k]):
                    sumErrorOver += diff**3
                    out[i][j][k][0] = 255
                    out[i][j][k][1] = -255*(diff/diffMax)+255
                    out[i][j][k][2] = -255*(diff/diffMax)+255
                    
                else:
                    out[i][j][k][0] = -255*(diff/diffMax)+255
                    out[i][j][k][1] = -255*(diff/diffMax)+255
                    out[i][j][k][2] = 255
                sumError += diff**3
    return out, sumError, sumErrorOver