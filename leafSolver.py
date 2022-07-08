# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:15:32 2019

@author: NickolasT
"""
from numba import jit
import numpy as np
from leaf import *
import translator as tr
import time
import grapher as gr
import matplotlib.pyplot as plt

# This file is used for the solver object, which stores general parameters
# for the solving algorithms
# It contains the current status of the virtual MLC--ie give the expected fluence
# and functions to execute particular solving algorithms that position the leaves

# This class could not be jitclass'd (due to use of leaf objects)
# so pertinant functions have been numba-ized at the bottom
class leafSolver:
    nLvls = 2                           # Number of layers to use
    xMin= -4.875                        # Minimum x-coord
    xMax = 4.875                        # Maximum x-coord
    xStep = 0.25                        # x-domain descretization step size
    step = 1E-4                         # Step size for searching values along x-domain
    maxWeight = 5                       # Maximum allow weight
    weightStep = 1E-2                   # Step size for search weights
    
    # Fraction of attenuation for layers. Indexing [n][i] where
    # n is the number of layers, and i is the layer (i=0 is the most attenuation)
    Tms = np.array([[0.02,0,0,0,0],[0.039,0.519,0,0,0],[0.062,0.446,0.723,0,0],[0.087,0.407,0.677,0.839,0],[0.112,0.378,0.641,0.812,0.907],[0.13256,0.35,0.62,0.8,0.89,0.93]],dtype=object)
    TmsPerm = []
    perms = []
    N = int((xMax-xMin)/xStep)          # Number of x descretizations
    leafLeft = []                       # List of leaf objects for left side
    leafRight = []                      # List of leaf objects for right side
    weight = []
    x = []                              # List of x positions
    nRows = 20                          #20 is 0.5  leaves.  10 is 1.0cm leaves
    
    def __init__(self, params):
        self.nLvls = params[0]
        self.xMin = params[1]
        self.xMax = params[2]
        self.xStep = params[3]
        self.step = params[4]
        self.maxWeight = params[5]
        self.weightStep = params[6]
            
        self.x = np.arange(self.xMin,self.xMax+self.xStep,self.xStep)
        self.N = len(self.x)
        
    def setup(self, rows):
        self.nRows = rows
        self.leafLeft = []
        self.leafRight = []
        # Instantiate the leaf objects
        for j in range(rows):
            temp = []
            for i in range(self.nLvls):
                temp.append(leaf(self.xMax,True,self.Tms[self.nLvls-1][i],self.xMin,self.xMax,i,j))
            self.leafLeft.append(temp)
        for j in range(rows):
            temp = []
            for i in range(self.nLvls):
                temp.append(leaf(self.xMin,False,self.Tms[self.nLvls-1][i],self.xMin,self.xMax,i,j))
            self.leafRight.append(temp)  
        self.nLvls = len(self.leafRight[0])
        
        # Calculates transmissions of all permutations of leaf coverage
        # Indexing is of binary 0=000, 1=001, 2=010, etc. and a binary 1 means
        # that leaf is closed
        # Where the LSB is the thinnest leaf, and MSB is the thickest Leaf
        self.TmsPerm = []
        for i in range(2**self.nLvls):
            temp = 1
            for j in range(self.nLvls-1,-1,-1):
                if(i%2 == 1):
                    temp *= self.Tms[self.nLvls-1][j]
                i = i//2
            self.TmsPerm.append(temp)
            
    def getFluence(self, weight):
        # Calculate the 2D array of the resulting fluence of the leaves
        # using a weight
        out = np.ones((len(self.leafLeft), self.N))
        for i in range(len(self.leafLeft)):
            out[i] *= self.getFluenceRow(weight,i)
        return out
    def getFluenceRow(self,weight,row):
        # Calculate the 1D array of the resulting fluence of the leaves
        # for a particular row using a weight
        out = np.ones(self.N) * weight
        for j in range(self.nLvls):             
            out *= self.leafLeft[row][j].mapTrans(self.N)
            out *= self.leafRight[row][j].mapTrans(self.N)
        return out
    
    def getLeaves(self):
        # Returns a list of the leaf positions rounded to 3 decimals
        left = []
        right = []
        for i in range(len(self.leafLeft)):
            temp = [] 
            for j in range(len(self.leafLeft[0])):
                temp.append(round(self.leafLeft[i][j].pos,3))
            left.append(temp)
        for i in range(len(self.leafRight)):
            temp = []
            for j in range(len(self.leafRight[0])):
                temp.append(round(self.leafRight[i][j].pos,3))
            right.append(temp)
        return (left, right)
    def getLeafTrans(self):
        out = []
        for L in self.leafLeft[0]:
            out.append(L.trans)
        return np.array(out)
    def setLeaves(self, left, right):
        # Assigns positions of leaves by 2D lists
        if not (len(left)==len(self.leafLeft) and len(right)==len(self.leafRight)
                and len(left[0])==len(self.leafLeft[0])) and len(right[0])==len(self.leafRight[0]):
            raise
        for i in range(len(left)):
            for j in range(len(left[0])):
                self.leafLeft[i][j].pos = left[i][j]
        for i in range(len(right)):
            for j in range(len(right[0])):
                self.leafRight[i][j].pos = right[i][j]
        return
#=================================================================
#*****************************************************************
#                       SOLVING METHODS
#*****************************************************************
#=================================================================
    def solveNone(self,inFluence):
        # Solve method that does not affect leaf positions or weight
        # Used for imported leaf positions
        inFluence = tr.interpolateData(inFluence,len(self.x))
        initTime = time.time()
        # Returns with the same output pattern as the other solve methods
        return (self.getFluence(self.weight), inFluence, self.x, time.time() - initTime)
    
    
    
    def solveExtend2PeakAndTrough(self,inFluence,troughThresh,peakFocus):
        # Algorithm to solve fluence pattern
        # Iterates through different weights to find lowest error
        # For each row, determine two peaks, left and right. Determine
        # permutation of leaves to attenuate to the peaks.
        # Position the other leaves to the outer edge of the respective peaks
        # extend permutated leaves to the opposite peak, assuming overtravel.
        # peakFocus is a boolean of which peak to focus on (True=upperPeak | False=lowerPeak)
        TROUGH_THRESH_FACTOR = troughThresh   # The most amount of trough conservation, factor relative
                                    # to the lower peak
        WEIGHT_LIMIT_FACTOR = 10      # Limits weight domain by [maxFluence/2,maxFluence*2]
        
        initTime = time.time()
        print("==============================")
        if peakFocus:
            print("Solving leaf positions with Extend2PeakAndTrough (Focus Upper Peak)")
        else:
            print("Solving leaf positions with Extend2PeakAndTrough (Focus Lower Peak)")
        print("------------------------------")
        self.setup(len(inFluence))                    # Instantiates all leaf objects
        inFluence = tr.interpolateData(inFluence,len(self.x))
        
        bestError = 10000        # Keeps track of best error
        
        # Iterate weights within range around global max weight and solve
        weights = np.arange(np.amax(inFluence)/WEIGHT_LIMIT_FACTOR,min([self.maxWeight,np.amax(inFluence)*WEIGHT_LIMIT_FACTOR]),self.weightStep)
        
        for w,weight in enumerate(weights):
            # Print status every 5% progress
            #if w % (len(weights)//10)==-0:
            #    print("Weight: %.4f / %.4f" % (weight, max(weights)))
                
            self.weight = weight
            
            # *******SOLVE THE LEAVES****** 
            # ---NUMBA FUNCTION SEE BELOW---
            L,R = NUMBA2PeakAndTrough(weight,self.x,self.nLvls,self.step,self.getLeafTrans(), inFluence,peakFocus,TROUGH_THRESH_FACTOR)             
            self.setLeaves(L,R)  
                    
            # Calculate the total error for the solution of this weight
            error = gr.diffMap(self.getFluence(self.weight), inFluence)[1]
            # Remember the best weight iteration
            if (bestError > error):
                bestError = error
                bestPos = self.getLeaves()
                bestWeight = self.weight
        # After iterating the weights, restore solver to have the best iteration
        for i in range(len(bestPos[0])):
            for j in range(len(bestPos[0][0])):
                self.leafLeft[i][j].pos = bestPos[0][i][j]
                self.leafRight[i][j].pos = bestPos[1][i][j]
        self.weight = bestWeight
        # Return the resulting fluence from the leaves
        # the inputted fluence--having been interpolated to the same x domain,
        # the x domain used, and the time elapsed
        return (self.getFluence(self.weight), inFluence, self.x, time.time() - initTime)
    
    
    def solveExtend2PeakWeightLimit(self,inFluence):
        # Algorithm to solve fluence pattern
        # Iterates through different weights to find lowest error
        # For each row, determine two peaks, left and right. Determine
        # permutation of leaves to attenuate to the peaks.
        # Position the other leaves to the outer edge of the respective peaks
        # extend permutated leaves to the opposite peak, assuming overtravel.
        WEIGHT_LIMIT_FACTOR = 10      # Limits weight domain by [maxFluence/2,maxFluence*2]
        
        initTime = time.time()

        self.setup(len(inFluence))                    # Instantiates all leaf objects
        inFluence = tr.interpolateData(inFluence,len(self.x))
        
        bestError = 10000        # Keeps track of best error
        
        # Iterate weights within range around global max weight and solve
        weights = np.arange(np.amax(inFluence)/WEIGHT_LIMIT_FACTOR,min([self.maxWeight,np.amax(inFluence)*WEIGHT_LIMIT_FACTOR]),self.weightStep)
                
        for w,weight in enumerate(weights):
            # Print status every 5% progress
            #if w % (len(weights)//10)==-0:
            #    print("Weight: %.4f / %.4f" % (weight, max(weights)))
                
            self.weight = weight
            # *******SOLVE THE LEAVES******
            # ---NUMBA FUNCTION SEE BELOW---
            L,R = NUMBA2PeakWeightLimit(weight,self.x,self.nLvls,self.step,self.getLeafTrans(), inFluence)             
            self.setLeaves(L,R)
            # Calculate the total error for the solution of this weight
            error = gr.diffMap(self.getFluence(self.weight), inFluence)[1]
            # Remember the best weight iteration
            if (bestError > error):
                bestError = error
                bestPos = self.getLeaves()
                bestWeight = self.weight
        # After iterating the weights, restore solver to have the best iteration
        for i in range(len(bestPos[0])):
            for j in range(len(bestPos[0][0])):
                self.leafLeft[i][j].pos = bestPos[0][i][j]
                self.leafRight[i][j].pos = bestPos[1][i][j]
        self.weight = bestWeight
        print("Best weight: "+str(bestWeight)[0:4])
        # Return the resulting fluence from the leaves
        # the inputted fluence--having been interpolated to the same x domain,
        # the x domain used, and the time elapsed
        return (self.getFluence(self.weight), inFluence, self.x, time.time() - initTime)
    
        
    
# ============================
# Numba-ized member functions
# =============================   
# leafSolver class was unable to be jitclass'd while utilizing the leaf object, so
# pertinant functions have been NUMBA'd below and called by the leafSolver accordingly
# Numba-izing the functions generally involve circumventing the use of leaf objects and
# handling only their specific indices and values
@jit(nopython=True)
def NUMBA2PeakAndTrough(weight,x,nLvls,step,leafTrans,inFluence,peakFocus,TROUGH_THRESH_FACTOR):
    N = len(x)
    leafLeftPos = np.zeros((len(inFluence),nLvls))
    leafRightPos = np.zeros((len(inFluence),nLvls))
    for r in range(len(inFluence)):
        # Find all maxima
        left = True
        maxFluence = 0
        maxL = 0
        maxR = 0
        maxXL = 0
        maxXR = 0
        for i in range(1,N-1):
            if (left and (inFluence[r][i-1] < inFluence[r][i] and inFluence[r][i+1] <= inFluence[r][i])):
                maxL = inFluence[r][i]
                maxXL = x[i]
                left = False
            if (inFluence[r][i-1] < inFluence[r][i] and inFluence[r][i+1] <= inFluence[r][i]):
                if inFluence[r][i] > maxFluence:
                    maxFluence = inFluence[r][i]
                    maxX = x[i]
                maxR = inFluence[r][i]
                maxXR = x[i]
        # Attempt to distinguish the major two peaks
        # ***THIS ALGORITHM MAY NEED TO BE IMPROVED***
        if (maxFluence != maxL and maxFluence != maxR):
            maxR = maxFluence
            maxXR = maxX
            
        # Initialize leaf positions at the peaks
        for i in range(nLvls):
            leafLeftPos[r][i] = maxXL
            leafRightPos[r][i] = maxXR
        # ****************************
        # Determine the permutation of leaves to attenuate to the
        # respective peak
        # ------------LEFT-----------
        for i in range(2**nLvls):
            prodL = weight
            fitL = np.ones(nLvls)*-1  # List of leaf indices to shape outer edge left peak
            extL = np.ones(nLvls)*-1  # List of leaf indices to shape inner edge right peak
            kFitL = 0
            kExtL = 0
            for j in range(nLvls-1,-1,-1):
                if (i%2 == 1):      # Use binary to iterate perms
                    prodL *= leafTrans[j]
                    extL[kExtL] = j
                    kExtL += 1
                else:
                    fitL[kFitL] = j
                    kFitL += 1
                i = i//2
            if (prodL <= maxL):
                break
                        
        # --------------RIGHT-------------
        for i in range(2**nLvls):
            prodR = weight
            fitR = np.ones(nLvls)*-1       # List of leaf indices to shape outer edge right peak
            extR = np.ones(nLvls)*-1       # List of leaf indices to shape inner edge left peak
            kFitR = 0
            kExtR = 0
            for j in range(nLvls-1,-1,-1):
                if (i%2 == 1):          # Use binary to iterate perms
                    prodR *= leafTrans[j]
                    extR[kExtR] = j
                    kExtR += 1
                else:
                    fitR[kFitR] = j
                    kFitR += 1
                i = i//2
            if (prodR <= maxR):
                break
        
        # ***Maintain trough by using leaves from unfocused peak to make trough***
        # It will try to maintain trough until fluence reaches the minimum between
        # the 2 peaks, OR fluence reaches a certain fraction (TROUGH_THRESH_FACTOR)
        # of the unfocused peak
        # TROUGH_THRESH_FACTOR is used to prevent over compensation to maintain
        # trough if the trough is very small relative to peaks
        # peakFocus : True=upperPeak | False=lowerPeak
        
        # Find minimum between the two peaks
        lowY = maxL
        for i in range(0,N):
            if x[i] > maxXR:
                break
            if (x[i] > maxXL) and (inFluence[r][i] < lowY):
                lowY = inFluence[r][i]
        # Assume fluence shielded by all leaves extended in between the two
        # peaks from initial leaf sorting
        prod = weight
        for i in range(len(extL)):
            if extL[i] >= 0:
                prod *= leafTrans[int(extL[i])]
        for i in range(len(extR)):
            if extR[i] >= 0:
                prod *= leafTrans[int(extR[i])]
        # Decide which leaves from unfocused peak to use instead for focused peak
        if peakFocus and (maxL > maxR) and len(fitR) > 0:       # Use upper peak
            while prod > lowY and prod > maxR * TROUGH_THRESH_FACTOR:
                j = 0
                # Find unimportant leaf to transfer
                while j < len(extR)-1 and leafTrans[int(extR[j])] > leafTrans[int(fitR[0])]:
                    j += 1
                # Transfer and insert the leaf into the extend group
                if j != nLvls:
                    for i in range(nLvls-1,j,-1):
                        extR[i] = extR[i-1]
                extR[j] = fitR[0]
                # Pop off the leaf from previous group
                fitR[0:-1] = fitR[1:]
                if fitR[0] < 0:
                    break
                prod *= leafTrans[int(fitR[0])]
        elif peakFocus and (maxR > maxL) and len(fitL) > 0 :    # Use upper peak
            while prod > lowY and prod > maxL * TROUGH_THRESH_FACTOR:
                j = 0
                # Find unimportant leaf to transfer
                while j < len(extL)-1 and leafTrans[int(extL[j])] > leafTrans[int(fitL[0])]:
                    j += 1 
                # Transfer and insert the leaf into the extend group
                if j != nLvls:
                    for i in range(nLvls-1,j,-1):
                        extL[i] = extL[i-1]
                extL[j] = fitL[0]
                # Pop off the leaf from previous group
                fitL[0:-1] = fitL[1:]
                if fitL[0] < 0:
                    break
                prod *= leafTrans[int(fitL[0])]
        elif not peakFocus and (maxL < maxR) and len(fitR) > 0:       # Use lower peak
            while prod > lowY and prod > maxR * TROUGH_THRESH_FACTOR:
                j = 0
                # Find unimportant leaf to transfer
                while j < len(extR)-1 and leafTrans[int(extR[j])] > leafTrans[int(fitR[0])]:
                    j += 1
                # Transfer and insert the leaf into the extend group
                if j != nLvls:
                    for i in range(nLvls-1,j,-1):
                        extR[i] = extR[i-1]
                extR[j] = fitR[0]
                # Pop off the leaf from previous group
                fitR[0:-1] = fitR[1:]
                
                if fitR[0] < 0:
                    break
                prod *= leafTrans[int(fitR[0])]
        elif not peakFocus and (maxR < maxL) and len(fitL) > 0 :    # Use lower peak
            while prod > lowY and prod > maxL * TROUGH_THRESH_FACTOR:
                j = 0
                # Find unimportant leaf to transfer
                while j < len(extL)-1 and leafTrans[int(extL[j])] > leafTrans[int(fitL[0])]:
                    j += 1
                # Transfer and insert the leaf into the extend group
                if j != nLvls:
                    for i in range(nLvls-1,j,-1):
                        extL[i] = extL[i-1]
                extL[j] = fitL[0]
                # Pop off the leaf from previous group
                fitL[0:-1] = fitL[1:]
                
                if fitL[0] < 0:
                    break
                prod *= leafTrans[int(fitL[0])]
        # ***Fit leaves of fitL to outer edge of left peak***
        xTemp = maxXL
        yTemp = 1
        prodFit,prodExt = prodL,prodL
        for i in range(len(fitL)):
            if fitL[i] < 0:
                continue
            targ = (leafTrans[int(fitL[i])]*prodFit + prodFit) / 2
            while yTemp >= targ:
                xTemp -= step
                yTemp = np.interp(xTemp,x,inFluence[r])
            leafLeftPos[r][int(fitL[i])] = xTemp
            prodFit *= leafTrans[int(fitL[i])]
        # Fit leaves of fitR to outer edge of right peak
        xTemp = maxXR
        yTemp = 1
        prodFit = prodR
        for i in range(len(fitR)):
            if fitR[i] < 0:
                continue
            targ = (leafTrans[int(fitR[i])]*prodFit + prodFit) / 2
            while yTemp >= targ:
                xTemp += step
                yTemp = np.interp(xTemp,x,inFluence[r])
            leafRightPos[r][int(fitR[i])]= xTemp
            prodFit *= leafTrans[int(fitR[i])]
        
        # *************************************
        # Extend the rest of the leaves until they meet the inner edge
        # of the opposite peak or the meet the opposite peak itself
        # ---------LEFT-------------
        prodExt = prodR
        
        xTemp = maxXR
        yTemp,yOld = 1,1
        for i in range(len(extL)):
            if extL [i] < 0:
                continue
            targ = (leafTrans[int(extL[i])]*prodExt + prodExt) / 2
            while yTemp >= targ or yOld <= yTemp:
                if (xTemp <= maxXL):
                    if (maxL > maxR):
                        xTemp = maxXR
                    break
                yOld = yTemp
                xTemp -= step
                yTemp = np.interp(xTemp, x, inFluence[r])

            leafLeftPos[r][int(extL[i])] = xTemp
            prodExt *= leafTrans[int(extL[i])]
        #------------RIGHT----------
        prodExt = prodL
        xTemp = maxXL
        yTemp, yOld = 1,1
        for i in range(len(extR)):
            if extR[i] < 0:
                continue
            targ = (leafTrans[int(extR[i])]*prodExt + prodExt) / 2
            while yTemp >= targ or yOld <= yTemp:
                if (xTemp >= maxXR):
                    if (maxR > maxL):
                        xTemp = maxXL
                    break
                yOld = yTemp
                xTemp += step
                yTemp = np.interp(xTemp, x, inFluence[r])
            leafRightPos[r][int(extR[i])] = xTemp
            prodExt *= leafTrans[int(extR[i])]
            
    return leafLeftPos, leafRightPos
        
@jit(nopython=True)
def NUMBA2PeakWeightLimit(weight,x,nLvls,step,leafTrans,inFluence):
    N = len(x)
    leafLeftPos = np.zeros((len(inFluence),nLvls))
    leafRightPos = np.zeros((len(inFluence),nLvls))
    for r in range(len(inFluence)):
        # Find all maxima
        left = True
        maxFluence = 0
        maxL = 0
        maxR = 0
        maxXL = 0
        maxXR = 0
        for i in range(1,N-1):
            if left and (inFluence[r][i-1] < inFluence[r][i] and inFluence[r][i+1] <= inFluence[r][i]):
                maxL = inFluence[r][i]
                maxXL = x[i]
                left = False
            if (inFluence[r][i-1] < inFluence[r][i] and inFluence[r][i+1] <= inFluence[r][i]):
                if inFluence[r][i] > maxFluence:
                    maxFluence = inFluence[r][i]
                    maxX = x[i]
                maxR = inFluence[r][i]
                maxXR = x[i]
        # Attempt to distinguish the major two peaks
        # ***THIS ALGORITHM MAY NEED TO BE IMPROVED***
        if (maxFluence != maxL and maxFluence != maxR):
            maxR = maxFluence
            maxXR = maxX
            
            
        # Initialize leaf positions at the peaks
        for i in range(nLvls):
            leafLeftPos[r][i] = maxXL
            leafRightPos[r][i] = maxXR
            
        # ****************************
        # Determine the permutation of leaves to attenuate to the
        # respective peak
        # ------------LEFT-----------
        for i in range(2**nLvls):
            prodL = weight
            fitL = np.ones(nLvls)*-1  # List of leaf indices to shape outer edge left peak
            extL = np.ones(nLvls)*-1  # List of leaf indices to shape inner edge right peak
            kFit = 0
            kExt = 0
            for j in range(nLvls-1,-1,-1):
                if (i%2 == 1):      # Use binary to iterate perms
                    prodL *= leafTrans[j]
                    extL[kExt] = j
                    kExt += 1
                else:
                    fitL[kFit] = j
                    kFit += 1
                i = i//2
            if (prodL <= maxL):
                break
                        
        # --------------RIGHT-------------
        for i in range(2**nLvls):
            prodR = weight
            fitR = np.ones(nLvls)*-1       # List of leaf indices to shape outer edge right peak
            extR = np.ones(nLvls)*-1       # List of leaf indices to shape inner edge left peak
            kFit = 0
            kExt = 0
            for j in range(nLvls-1,-1,-1):
                if (i%2 == 1):          # Use binary to iterate perms
                    prodR *= leafTrans[j]
                    extR[kExt] = j
                    kExt += 1
                else:
                    fitR[kFit] = j
                    kFit += 1
                i = i//2
            if (prodR <= maxR):
                break
        
        # Fit leaves of fitL to outer edge of left peak
        xTemp = maxXL
        yTemp = 1
        prodFit,prodExt = prodL,prodL
        for i in range(len(fitL)):
            if fitL[i] < 0:
                continue
            targ = (leafTrans[int(fitL[i])]*prodFit + prodFit) / 2
            while yTemp >= targ:
                xTemp -= step
                yTemp = np.interp(xTemp,x,inFluence[r])
            leafLeftPos[r][int(fitL[i])] = xTemp
            prodFit *= leafTrans[int(fitL[i])]
        # Fit leaves of fitR to outer edge of right peak
        xTemp = maxXR
        yTemp = 1
        prodFit = prodR
        for i in range(len(fitR)):
            if fitR[i] < 0:
                continue
            targ = (leafTrans[int(fitR[i])]*prodFit + prodFit) / 2
            while yTemp >= targ:
                xTemp += step
                yTemp = np.interp(xTemp,x,inFluence[r])
            leafRightPos[r][int(fitR[i])]= xTemp
            prodFit *= leafTrans[int(fitR[i])]
        
        # *************************************
        # Extend the rest of the leaves until they meet the inner edge
        # of the opposite peak or the meet the opposite peak itself
        # ---------LEFT-------------
        prodExt = prodR
        
        xTemp = maxXR
        yTemp,yOld = 1,1
        for i in range(len(extL)):
            if extL [i] < 0:
                continue
            targ = (leafTrans[int(extL[i])]*prodExt + prodExt) / 2
            while yTemp >= targ or yOld <= yTemp:
                if (xTemp <= maxXL):
                    if (maxL > maxR):
                        xTemp = maxXR
                    break
                yOld = yTemp
                xTemp -= step
                yTemp = np.interp(xTemp, x, inFluence[r])

            leafLeftPos[r][int(extL[i])] = xTemp
            prodExt *= leafTrans[int(extL[i])]
        #------------RIGHT----------
        prodExt = prodL
        xTemp = maxXL
        yTemp, yOld = 1,1
        for i in range(len(extR)):
            if extR[i] < 0:
                continue
            targ = (leafTrans[int(extR[i])]*prodExt + prodExt) / 2
            while yTemp >= targ or yOld <= yTemp:
                if (xTemp >= maxXR):
                    if (maxR > maxL):
                        xTemp = maxXL
                    break
                yOld = yTemp
                xTemp += step
                yTemp = np.interp(xTemp, x, inFluence[r])
            leafRightPos[r][int(extR[i])] = xTemp
            prodExt *= leafTrans[int(extR[i])]
            
    return leafLeftPos, leafRightPos
        