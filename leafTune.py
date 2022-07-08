# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:50:17 2019

@author: anmorrow
"""
import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import radon
from copy import deepcopy
import mapper as mp
from leafSolver import *

def leaftune(params, solverParams, theta, origSolverArray,origWeights, volume0, initialize):
    #start up
    #theta = np.linspace(-179, 179,178, endpoint=False)
    #theta = [0,40,80,120,160,200,240,280,320]
    nRows = len(origSolverArray[0].leafLeft)
    nLvls= params[0]
    
    newSolverArray = solverCopy(origSolverArray, params, nRows, nLvls)
    
    #sinoIdealAll = np.zeros((nRows,178,40))#leaf row, angle, x position 
    sinoIdealAll = np.zeros((nRows,len(theta),40))#leaf row, angle, x position 
    
    for i in range(len(volume0)):
        sinoIdealAll[i] = np.transpose(radon(volume0[i], theta=theta))
    newWeights = deepcopy(origWeights)
    newfluencesArray = mp.getSolverArrayFluencesWithWeights(newSolverArray,newWeights)
    newVolume1 = mp.mapFluences(newfluencesArray,solverParams)
    subScore = mp.diffVolume(newVolume1,volume0)[1]

    if not(initialize):

        
        print(subScore)
        
        #optimize weights
        oldfluencesArray = deepcopy(newfluencesArray)
        oldSolverArray = solverCopy(newSolverArray, params, nRows, nLvls)
        oldsubScore = subScore
        origScore = subScore
        weightStepArray=[1.01,0.99]
        turn = 0
        
        while True:
            origScore = subScore
            changearr=len(theta)*[0]
            print(turn)
            for ang in range(len(theta)):
                thisStepArray = deepcopy(weightStepArray)
                while len(thisStepArray)>0:
                    oldfluencesArray = deepcopy(newfluencesArray)
                    oldSolverArray = solverCopy(newSolverArray, params, nRows, nLvls)
                    oldsubScore = subScore
                    newSolverArray[ang].weight*=thisStepArray[0]
                    newWeights[ang] =newSolverArray[ang].weight
                    if newSolverArray[ang]==0:
                        newSolverArray[ang].weight=thisStepArray[0]
                        newWeights[ang] =newSolverArray[ang].weight
                    newfluencesArray = mp.getSolverArrayFluencesWithWeights(newSolverArray,newWeights)
                    newVolume1 = mp.mapFluences(newfluencesArray,solverParams)
                    subScore = mp.diffVolume(newVolume1,volume0)[1]
                    
                    if subScore>=oldsubScore:
                        newfluencesArray=deepcopy(oldfluencesArray)
                        newSolverArray = solverCopy(oldSolverArray, params, nRows, nLvls)
                        subScore = oldsubScore
                        print("Fld "+str(ang) + " weight change DID NOT WORK. Score: "+str(subScore)[:7])
                        del thisStepArray[0]                    
                    else:
                        print("Fld "+str(ang) + " "+str(thisStepArray[0])+" weight change WORKED. Score: " +str(subScore)[:7]) 
                        changearr[ang]=1
            turn+=1
            if subScore == origScore or sum(changearr)==0 or sum(changearr)==1:#if there is no change in score, or only 1 or 0 fields change, stop
                break
            
        #optimize leaf positions    
        #neglect slices that have no fluence
        rowsToSearch=[]
        for row in range(nRows):
            if sum(sum(volume0[row]))!=0:
                rowsToSearch.append(row)
            else:
                for ang in range(len(theta)):
                    for lvl in range(nLvls):
                        newSolverArray[ang].leafRight[row][lvl].pos=0
                        newSolverArray[ang].leafLeft[row][lvl].pos=0
                        oldSolverArray[ang].leafRight[row][lvl].pos=0
                        oldSolverArray[ang].leafLeft[row][lvl].pos=0
    

        for ang in range(len(theta)):
            #future improvement - calc volume fluence for all but the current angle being optimized.  add it back in with mapfluences
            for lvl in range(nLvls):
                #right side   
                for row in rowsToSearch:
                    leafPos = newSolverArray[ang].leafRight[row][lvl].pos
                    pos=leafPos
                    xPos = round((leafPos+4.875)/0.25)
                    sinoErr=np.transpose(radon(newVolume1[row], theta=[theta[ang]]))[0]
                    sinoErrVal = sinoErr[xPos]-sinoIdealAll[row, ang,xPos]
                    while xPos !=0 and xPos != 40:
                        if sinoErrVal > 0:
                            xPos = xPos-1
                            if xPos<0:
                                break
                            if (sinoErr[xPos] ==0 or sinoErr[xPos]>0):
                                break
                            else:
                                pos = xPos/4.0 - 4.875
                                newSolverArray[ang].leafRight[row][lvl].pos=pos
                        elif sinoErrVal < 0:
                            xPos = xPos+1
                            if xPos >39:
                                break
                            if (sinoErr[xPos] ==0 or sinoErr[xPos]<0):
                                break
                            else:
                                pos = xPos/4.0 - 4.875
                                newSolverArray[ang].leafRight[row][lvl].pos=pos
                        elif sinoErrVal == 0:
                            break
                    if leafPos != pos:
                        newfluencesArray = mp.getSolverArrayFluencesWithWeights(newSolverArray,newWeights)
                        newVolume1 = mp.mapFluences(newfluencesArray,solverParams)
                        subScore = mp.diffVolume(newVolume1,volume0)[1]
                        if subScore>=oldsubScore:
                            newfluencesArray=deepcopy(oldfluencesArray)
                            newSolverArray = solverCopy(oldSolverArray, params, nRows, nLvls)
                            subScore = oldsubScore
                            print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Rt leaf change DID NOT WORK. Score: "+str(subScore)[:7])
                            
                        else:
                            oldfluencesArray = deepcopy(newfluencesArray)
                            oldSolverArray = solverCopy(newSolverArray, params, nRows, nLvls)
                            oldsubScore = subScore
                            print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Rt leaf change WORKED. Score: "+str(subScore)[:7]) 
                    elif leafPos == pos:
                        newfluencesArray=deepcopy(oldfluencesArray)
                        newSolverArray = solverCopy(oldSolverArray, params, nRows, nLvls)
                        subScore = oldsubScore
                        print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Rt leaf NOT MOVED. Score: "+str(subScore)[:7])
                
                #left side
                
                for row in rowsToSearch:     
                    leafPos = newSolverArray[ang].leafRight[row][lvl].pos
                    pos=leafPos
                    xPos = round((leafPos+4.875)/0.25)
                    sinoErr=np.transpose(radon(newVolume1[row], theta=[theta[ang]]))[0]
                    sinoErrVal = sinoErr[xPos]-sinoIdealAll[row, ang,xPos]
                    while xPos !=0 and xPos != 40:
                        if sinoErrVal < 0:
                            xPos = xPos-1
                            if xPos < 0:
                                break
                            if (sinoErr[xPos] ==0 or sinoErr[xPos]>0):
                                break
                            else:
                                pos = xPos/4.0 - 4.875
                                newSolverArray[ang].leafLeft[row][lvl].pos=pos
                        elif sinoErrVal >0:
                            xPos = xPos+1
                            if xPos >39:
                                break
                            if (sinoErr[xPos] ==0 or sinoErr[xPos]<0):
                                break
                            else:
                                pos = xPos/4.0 - 4.875
                                newSolverArray[ang].leafLeft[row][lvl].pos=pos
                        elif sinoErrVal == 0:
                            break
                    if leafPos != pos:
                        
                        newfluencesArray = mp.getSolverArrayFluencesWithWeights(newSolverArray,newWeights)
                        newVolume1 = mp.mapFluences(newfluencesArray,solverParams)
                        subScore = mp.diffVolume(newVolume1,volume0)[1]
                        
                        if subScore>=oldsubScore:
                            newfluencesArray=deepcopy(oldfluencesArray)
                            newSolverArray = solverCopy(oldSolverArray, params, nRows, nLvls)
                            subScore = oldsubScore
                            print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Lt leaf change DID NOT WORK. Score: "+str(subScore)[:7])
                            
                        else:
                            oldfluencesArray = deepcopy(newfluencesArray)
                            oldSolverArray = solverCopy(newSolverArray, params, nRows, nLvls)
                            oldsubScore = subScore
                            print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Lt leaf change WORKED. Score: "+str(subScore)[:7]) 
                    elif leafPos == pos:
                        newfluencesArray=deepcopy(oldfluencesArray)
                        newSolverArray = solverCopy(oldSolverArray, params, nRows, nLvls)
                        subScore = oldsubScore
                        print("Fld "+str(ang) + ", lvl "+str(lvl) + ", row " +str(row) + " Lt leaf NOT MOVED. Score: "+str(subScore)[:7])
            
        # plt.imshow(sinoIdealAll[10])
        # plt.imshow(sinoErr)
        # plt.plot(sinoErr[10])
        # gr.drawLeaves(newSolverArray[10],10)

# Draw the final resulting volume
#diffMax = np.amax(np.abs(newVolume1-volume0))
#fluenceMax = np.amax(np.vstack([newVolume1,volume0]))
#for i in range(nRows):
#        fig,_,_,_ = gr.mapFluenceREL(x,volume0[i],newVolume1[i],(0,np.amax(newVolume1[i])),diffMax,fluenceMax)
        #fig.savefig("row"+str(i)+"_"+str(nLvls)+"lvls[FINAL].png")
    return newSolverArray, newWeights, newVolume1, subScore
        
     
def blankSolver(solverToCopy, params, nRows):
    test=[]
    for i in range(len(solverToCopy)):
        test.append(leafSolver(params))
        test[i].setup(nRows)
    return test  
    
def solverCopy(solverToCopy, tparams, tnRows, tnLvls): #gotta make new versions of stuff otherwise were pointing at original memory positions
    test=[]
    for tSlice in range(len(solverToCopy)):
        
        test.append(leafSolver(tparams))
        test[tSlice].setup(tnRows)
        test[tSlice].weight = solverToCopy[tSlice].weight
        for tLeaf in range(len(solverToCopy[0].leafLeft)):
            for tLayer in range(tnLvls):
                test[tSlice].leafLeft[tLeaf][tLayer].pos=solverToCopy[tSlice].leafLeft[tLeaf][tLayer].pos
                test[tSlice].leafRight[tLeaf][tLayer].pos=solverToCopy[tSlice].leafRight[tLeaf][tLayer].pos
    
    return test