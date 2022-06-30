# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 08:58:58 2019

@author: anmorrow
"""

import pydicom
import numpy as np
from os import listdir
from os.path import isfile, join

import re
from math import log

dp = pydicom.dcmread(r"C:\Users\anmorrow.BHCS\Desktop\MLC_project\OLD\optimalFluence[2019_10_15]\RSTemplate.dcm")

mypath = r"C:\Users\anmorrow.BHCS\Desktop\MLC_project\optimalfluence[2019_10_15]\RS_PRS_D4_1LVL5mm\fluences"   #change this path per plan
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]             #get the list of fluences






for i in range(len(onlyfiles)):
    fluenceNum = onlyfiles[i][6:-16] #this says which field this one is.       
    fName = mypath + r"\Field " + str(fluenceNum) + ".optimal_fluence"  
    flFile = open(fName,"r")
    fl=flFile.readlines()
    fl=fl[9:49]  #just want the numerical data, 40x40
    flArray = []
    for j in range(len(fl)):
        row = np.asarray(re.split(r'\t+',fl[j])[:-1]).astype(np.float) #split with \t, get rit of \n at end of each line.  convert to a float array
        flArray=np.append(flArray,row)
    flFile.close() #got the data in the same format as the data is in the dicomRT file!  
    thickArray = [abs(-10*log(x)/0.5538) for x in flArray] #now we also have the compensator thickness array.  abs is to change -0.0 to 0.0
    #flArray = [x for x in flArray] #i dont know how stuff works so this is how i fix it
    flArray = [x/100 for x in flArray] # i think this is causing problems so i want to try this 
    
    dp.BeamSequence[int(fluenceNum)-1].CompensatorSequence[0].CompensatorThicknessData = thickArray #zero indexing correction
    dp.BeamSequence[int(fluenceNum)-1].CompensatorSequence[0].CompensatorTransmissionData = flArray #zero indexing correction
    
    #JUST FOR NOW:
    #dp.BeamSequence[i].ControlPointSequence[0].BeamLimitingDevicePositionSequence[0].LeafJawPositions[0] = "-40"
    #dp.BeamSequence[i].ControlPointSequence[0].BeamLimitingDevicePositionSequence[0].LeafJawPositions[1] = "40"
    
    #dp.BeamSequence[i].ControlPointSequence[0].BeamLimitingDevicePositionSequence[1].LeafJawPositions[0] = "-45"
    #dp.BeamSequence[i].ControlPointSequence[0].BeamLimitingDevicePositionSequence[1].LeafJawPositions[1] = "45"
    
#now export the file
dp.save_as(mypath+r"\newDicom.dcm")