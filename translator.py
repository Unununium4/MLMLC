# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:08:17 2019

@author: NickolasT
"""

import numpy as np
import pydicom
from pydicom.datadict import dictionary_VR
import mapper as mp
import random
from pydicom.dataset import Dataset

from os import listdir
from os.path import isfile, join
# This file parses a fluence file into the array of desired fluence values
# It includes some methods that can reshape or simplify the data
# It also exports details of solver object (leaf positions and weight)

def rewriteWeights(origfile, origfolder, newweights, fname):

    #first let's get the bounds of the RD array that we want to optimize with.  No need to optimize everywhere

    #get the mlmlc field doses
    count=0
    for f in listdir(origfolder):
        if isfile(join(origfolder, f)):
            if f[0:2]=="RD":
                tempd = pydicom.dcmread(origfolder +"\\"+f)
                scale = tempd.DoseGridScaling
                if count == 0:
                    dosetowrite = tempd.pixel_array * scale * newweights[count]
                else:
                    dosetowrite += tempd.pixel_array * scale* newweights[count]
                count += 1
            if f[0:2]=="RP":
                rp = pydicom.dcmread(origfolder +"\\"+f)
                rp.RTPlanLabel=fname
                rp[0x20,0xe].value = str(random.randint(0,1000000000000000000000000000000))#series instance id 
                rp.StudyInstanceUID = str(random.randint(0,1000000000000000000000000000000))
                rp.SOPInstanceUID = str(random.randint(0,1000000000000000000000000000000))
                #rp.FractionGroupSequence[0].ReferencedBeamSequence[b].BeamMeterset
                for b in range(len(rp.FractionGroupSequence[0].ReferencedBeamSequence)):
                    rp.FractionGroupSequence[0].ReferencedBeamSequence[b].BeamMeterset *= newweights[b]
                rp.save_as(fname+"plan.dcm")
    
    dosefile = pydicom.dcmread(origfile)
    scaling = dosefile.DoseGridScaling
    newarray = dosefile.pixel_array
    rows = len(newarray)
    columns = len(newarray[0])
    slices = len(newarray[0,0])


    for i in range(rows):
        for j in range(columns):
            for k in range (slices):
                newarray[i,j,k] =  np.uint32(dosetowrite[i,j,k] / scaling)
                
    newf = fname + "dose.dcm"
    
    print(np.max(newarray))
    dosefile.PixelData = newarray.tobytes()
    dosefile.SeriesInstanceID = str(random.randint(0,1000000000000000000000000000000))
    dosefile[0x20,0xe].value = str(random.randint(0,1000000000000000000000000000000))
    dosefile.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = str(random.randint(0,1000000000000000000000000000000))
    dosefile.StudyInstanceUID = str(random.randint(0,1000000000000000000000000000000))
    dosefile.SOPInstanceUID = str(random.randint(0,1000000000000000000000000000000))
    
    dosefile.SeriesDescription =newf
    dosefile.save_as(newf)
    
def convertToComp(plan, solverarray, fname):
    #change UIDs so we can reimport it
    plan.RTPlanLabel+="BadComp"
    plan[0x20,0xe].value = str(random.randint(0,1000000000000000000000000000000))#series instance id 
    plan.StudyInstanceUID = str(random.randint(0,1000000000000000000000000000000))
    plan.SOPInstanceUID = str(random.randint(0,1000000000000000000000000000000))
    
    mlmlcfluences=mp.getSolverArrayFluences(solverarray)
    
    cp1=Dataset()
    cp1.ControlPointIndex = "1"
    cp1.CumulativeMetersetWeight = "1.0"
    
    for b in range(len(plan.BeamSequence)):
            #compensator parts
        #start off by changing the dicom bits that already exist
        plan.BeamSequence[b].CompensatorSequence[0].MaterialID = 'Compensator'
        plan.BeamSequence[b].CompensatorSequence[0].SourceToCompensatorTrayDistance = 698.0
        plan.BeamSequence[b][0x300a, 0x00e0].value = 1 #num of compensators
        
        expandedArray = np.zeros((2*len(mlmlcfluences[b]),len(mlmlcfluences[b][0])))
        for i in range(len(mlmlcfluences[b])):
#           expandedArray[4*i]=fluenceArray[i]     #use the 4* for 1cm leaves
#           expandedArray[4*i+1]=fluenceArray[i]
#           expandedArray[4*i+2]=fluenceArray[i]
#           expandedArray[4*i+3]=fluenceArray[i]
            expandedArray[2*i]= mlmlcfluences[b][i]#use the 2* for the 0.5cm leaves
            expandedArray[2*i+1]= mlmlcfluences[b][i]
        fluences = expandedArray.flatten()
        
        for f in range(len(fluences)):#we are just gonna assume 0.5cm leaves so i have to double up the rows
            plan.BeamSequence[b].CompensatorSequence[0].CompensatorTransmissionData[f] = fluences[f]
        plan.BeamSequence[b].CompensatorSequence[0].CompensatorType = 'STANDARD'
        #add stuff that doesnt exist
        plan.BeamSequence[b].CompensatorSequence[0].CompensatorID = 'Compensator'
        plan.BeamSequence[b].CompensatorSequence[0].CompensatorDivergence = 'PRESENT'
        plan.BeamSequence[b].CompensatorSequence[0].CompensatorMountingPosition = 'SOURCE_SIDE'
        plan.BeamSequence[b].CompensatorSequence[0].add_new([0x300a,0x00ec], 'DS', plan.BeamSequence[b].CompensatorSequence[0].CompensatorTransmissionData)
        for f in range(len(fluences)):
            plan.BeamSequence[b].CompensatorSequence[0].CompensatorThicknessData[f] = abs(-10*np.log(fluences[f])/0.5538)
        #MLC parts
            #remove second part of compensator sequence
        del plan.BeamSequence[b].CompensatorSequence[1]
        
        #now get rid of the MLC control points after the first one
        for i in range(len(plan.BeamSequence[b].ControlPointSequence)-1):
            del plan.BeamSequence[b].ControlPointSequence[-1]

        #remove the mlc from that 1st cp that we left in, and a little something else
        del plan.BeamSequence[b].ControlPointSequence[0].BeamLimitingDevicePositionSequence[2]
        del plan.BeamSequence[b].ControlPointSequence[0].ReferencedDoseReferenceSequence
        del plan.BeamSequence[b].BeamLimitingDeviceSequence[2]
        #add back in a little bit that's needed into the cpsequence
        plan.BeamSequence[b].ControlPointSequence+=[cp1]
        #change number of cps
        plan.BeamSequence[b][0x300a, 0x0110].value = 2
    plan.save_as(fname)

def convertCPtoFluence(controlpoint, cpmu, xstep,xmin, xmax, ymin, ymax, dlg):
    #first figure out the size of the array we'll need
    nx = np.int16(np.round(abs(np.round((xmax+1-xmin)/xstep))))
    ny = np.int16(np.round(abs(np.round((ymax+1-ymin)/5))))
    array=np.zeros(nx, ny)
    firstyindex = 10+np.round((100+ymin)/5)#first 10 leaves are 10mm leaves.  we are pretending they dont exist except when we cant
    lasttyindex = 10+np.round(ymax/5)
    for y,leaf in enumerate(range(firstyindex,lasttyindex)):#120 leaves
        for x in range (nx) :
            xpos=xmin+x*xstep 
            if xpos > controlpoint[leaf] - dlg and xpos < controlpoint[leaf+60] + dlg:
                array[x,y]=cpmu
    return(array)

def importDicomFluence(beam):
    # Parses the fluence from the rt dicom file 
    index = []#index is going to be the gantry angle now
    sizeX = []
    sizeY = []
    spaceX = []
    spaceY = []
    originX = []
    originY = []
    data = []
    
    index = float(beam.ControlPointSequence[0].GantryAngle)
    sizeX = int(beam.CompensatorSequence[0].CompensatorColumns)
    sizeY = int(beam.CompensatorSequence[0].CompensatorRows)
    spaceX = float(beam.CompensatorSequence[0].CompensatorPixelSpacing[0])
    spaceY = float(beam.CompensatorSequence[0].CompensatorPixelSpacing[1])
    originX = float(beam.CompensatorSequence[0].CompensatorPosition[0])
    originY = float(beam.CompensatorSequence[0].CompensatorPosition[1])
    tempdata = beam.CompensatorSequence[0].CompensatorTransmissionData
    
    data = np.zeros((sizeY,sizeX))
    
        # unflatten and set the minimums to zero.  helps the algorithm along.
    minf = min(tempdata)
    for x in range(sizeX):
        for y in range(sizeY):
            if tempdata[y*sizeX + x] != minf:
                data[y,x] = float(tempdata[y*sizeX + x])

    return (index,sizeX,sizeY,spaceX,spaceY,originX,originY,data)


def importFluence(filename):
    # Parses the .optimal_fluence file 
    index = []
    sizeX = []
    sizeY = []
    spaceX = []
    spaceY = []
    originX = []
    originY = []
    data = []
    
    fd = open(filename,'r')
    for i,line in enumerate(fd.readlines()):
        if i == 0:
            index = line.split()[-3]
        elif i == 2:
            sizeX = int(line.split()[-1])
        elif i == 3:
            sizeY = int(line.split()[-1])
        elif i == 4:
            spaceX = float(line.split()[-1])
        elif i == 5:
            spaceY = float(line.split()[-1])
        elif i == 6:
            originX = float(line.split()[-1])
        elif i == 7:
            originY = float(line.split()[-1])
        elif i > 8:
            temp = []
            for j in line.split():
                temp.append(float(j))
            data.append(temp)
    if (len(data[-1])==0):
        data = data[:-1]
    fd.close()
    return (index,sizeX,sizeY,spaceX,spaceY,originX,originY,data)

def importFile(filename, delim, params):
    # Imports file containing leaf positions and returns solver    
    import leafSolver
    fd = open(filename,'r')
    lines = fd.readlines()
    
    mySolver = leafSolver.leafSolver(params)
    mySolver.nLvls = (len(lines[0].split(delim))-1) // 2
    mySolver.setup(len(lines)-1)
    for i,line in enumerate(lines):
        if i > 0:
            entries = line.split(delim)
            for j in range((len(entries)-1) // 2):
                mySolver.leafLeft[i-1][(len(entries)-1) // 2-1-j].pos = float(entries[1+j*2])
                mySolver.leafRight[i-1][(len(entries)-1) // 2-1-j].pos = float(entries[2+j*2])
        if i == 1:
            mySolver.weight = float(entries[-1])
    fd.close()
    return mySolver

def outputLeaves(filename, mySolver,delim):
    leafHeaders = ["Fluence Row", "Lt Leaf %d", "Rt Leaf %d", "weight"]
    
    fd = open(filename,'w')
    
    out = leafHeaders[0]
    for i in range(mySolver.nLvls):
        out += "\t" + (leafHeaders[1] % (i+1))
        out += "\t" + (leafHeaders[2] % (i+1))
    out += "\t" + leafHeaders[3]
    fd.write(out + "\n")
    
    for i in range(len(mySolver.leafLeft)):
        out = str(i+1)
        for j in range(len(mySolver.leafLeft[0])-1,-1,-1):
            out += "\t" + str(mySolver.leafLeft[i][j].pos)
            out += "\t" + str(mySolver.leafRight[i][j].pos)
        if i == 0:
            out += "\t" + str(mySolver.weight)
        fd.write(out + "\n")
    fd.close()
    
#ANM made this to export fluence files
def exportFluence(filename, field, fluenceArray):

    expandedArray = np.zeros((40,40))
    for i in range(len(fluenceArray)):
#        expandedArray[4*i]=fluenceArray[i]     #use the 4* for 1cm leaves
#        expandedArray[4*i+1]=fluenceArray[i]
#        expandedArray[4*i+2]=fluenceArray[i]
#        expandedArray[4*i+3]=fluenceArray[i]
        expandedArray[2*i]=fluenceArray[i]#use the 2* for the 0.5cm leaves
        expandedArray[2*i+1]=fluenceArray[i]
    
    fd = open(filename,'w')
    out = "#Field " + str(field) + " - fluence\n"
    out+="optimalfluence\n"
    out+="sizex\t40\n"
    out+="sizey\t40\n"
    out+="spacingx\t2.5\n"
    out+="spacingy\t2.5\n"
    out+="originx\t-48.75\n"
    out+="originy\t48.75\n"
    out+="values\n"
    #fd.write(out + "\n")
    
    for i in range(len(expandedArray)):
        for j in range(len(expandedArray[0])):
            out += str(expandedArray[i][j])+"\t"
        #fd.write(out + "\n")
        out+="\n"
    fd.write(out)
    fd.close()

    
def exportWeights(weights,normFluences ):
    fd = open("weights.txt",'w')
    out=""
    for w in range(len(weights)):
        out+=str(weights[w]) +"\n"
    fd.write(out)
    fd.close
    
    maxFluences = open("max_fluences.txt",'w')
    out=""
    for w in range(len(normFluences)):
        out+=str(np.max(normFluences[w])) +"\n"
    maxFluences.write(out)
    maxFluences.close
    
    
def normalizedKernel(data, n):
    # Simplifies the data by averaging every n-rows and then normalizes it
    mMax = np.amax(data)
    out = []
    for i in range(0,len(data),n):
        temp = []
        for k in range(len(data[0])):
            mSum = 0
            for j in range(n):
                mSum += data[i+j][k]
            temp.append(mSum / mMax / n)
        out.append(temp)
    
    return np.array(out)

def interpolateData(data,xN):
    # Increases the resolution of the data via interpolation
    # It will linearly interpolate until there are xN elements for the row
    xOldN = len(data[0])
    
    yN = len(data)                  # Temporary only interpolate in x-direction
    xNew = np.arange(0,xOldN,xOldN / xN)
    xOld = np.arange(0,xOldN+1)
    dataNew = np.ones((yN,xN))
    for j in range(yN):
        dataOld = data[j]
        dataOld = np.append(dataOld,data[j][-1])
        for i in range(xN):
            dataNew[j][i] = np.interp(xNew[i],xOld,dataOld)
    return dataNew
            