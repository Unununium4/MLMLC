# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:00:05 2019

@author: anmorrow
"""
import re
import numpy as np
from os import listdir
from os.path import isfile, join
import pydicom

def timetodeliver(plan, drate):
    ttd=0#time to deliver
    #assume 600MU/minute.  this is in MU/second - make this an input to allow for FFF rates
    lspeed = 25#assume leaf speed of 25mm per second projected at isocenter
    gspeed = 6#assume 6 degrees per second gantry rotation (360degree/minute)


    #get the list of MUs per control point.  
    #get the list of furthest that any leaf has to travel per control point
    #get the list of gantry angles that each control point has to cover
    #pick the longest time of the three per control point and add it up
    return ttd


#old way
# #change the two variables below per plan
# myPath = r"C:\Users\anmorrow.BHCS\Desktop\MLC_project\optimalfluence[2019_10_15]\RS_PRS_D4_6LVL5mm\leaves"
# mufile = open(r"C:\Users\anmorrow.BHCS\Desktop\MLC_project\optimalfluence[2019_10_15]\6Lvl5mmNW.txt","r")
# mus = mufile.readlines()

# #get num leaf files, number of leaf rows, and number of leaf levels to build an array
# leafSpeed = 2.5
# onlyFiles = [f for f in listdir(myPath) if isfile(join(myPath, f))]  #file list
# fName = myPath + "\\" +onlyFiles[0]
# leafFile = open(fName,"r")
# leaves = leafFile.readlines()
# leaves = leaves[1:] #first row is the header
# lastLeafPos= len(re.split(r'\t+',leaves[0]))-1 #first row has the weight of the field as the last column.  no other rows have this
# numFiles = len(onlyFiles)
# numRows = len(leaves)
# numLeaves = lastLeafPos - 1
# leafPosArray=np.zeros((numFiles,numRows,numLeaves))
# fileNameEnd = onlyFiles[0][-20:]# = _Nlvls[FINAL].leaves

# muTimeArray = np.zeros(numFiles) #MU/control point
# gantryTimeArray = np.zeros(numFiles)


# #start iterating through all.leaf files
# for i in range(numFiles):
#     fName = myPath + "\\" + str(i+1) + fileNameEnd#must read in files in their delivered order, not however listdir  orders them
#     leafFile = open(fName,"r")
#     leaves = leafFile.readlines()
#     leaves = leaves[1:]
#     #start iterating through all rows
#     for j in range(len(leaves)):
#         leafPosArray[i,j,:]=np.asarray(re.split(r'\t+',leaves[j])[1:lastLeafPos]).astype(np.float)
#     leafFile.close()
#     muTimeArray[i]=float(re.split(r'\t+',mus[i])[2][:-1])/10 #split by tabs, get rid of trailing /n
#     gantryTimeArray[i]=60.0/178.0#time to go between two control points
# #now that we have all of the leaf file data we need to find the time of travel for each leaf between each sequential control point
# #then find the max time of travel between each control point
# #then sum this over all control points
# totalTime = 0  #start out with total time = time to deliver the MU of the plan


# leafTimeArray=np.zeros((numFiles-1,numRows,numLeaves))

# for i in range(numFiles-1):
#     for j in range(numRows):
#         for k in range(numLeaves):
#             leafTimeArray[i,j,k]=abs(leafPosArray[i,j,k]-leafPosArray[i+1,j,k])/leafSpeed
# #cpTimes = zeros(numFiles-1)
# for i in range(numFiles-1):
#     totalTime += max([np.max(leafTimeArray[i]),gantryTimeArray[i],muTimeArray[i]])# pick the max time of gantry motion, leaf motion, or mu delivery.looks like we'll ignore that last little control for MU and gantry motion (<1s error likely)  
# print(totalTime)