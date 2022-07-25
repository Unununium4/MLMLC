# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 14:17:28 2022

@author: ANMORROW
"""

import pydicom
import numpy as np
from os import listdir
from os.path import isfile, join
from copy import deepcopy
import translator as tr
import time

#cost function
def cost(ideal, current):
    c = 100*np.sum(abs(ideal-current)/np.max(ideal))/(len(ideal)*len(ideal[0])*len(ideal[0][0]))
    return c #average % difference per pixel, relative to the maximum dose in the ideal set

def sumall(fielddoses, wts):
    s=fielddoses[0]*0
    for b in range(len(wts)):
        s+=fielddoses[b]*wts[b]
    return s

def IsPointInPolygon(contourxs,contourys, pointx, pointy):
    result = False
    j = len(contourxs) - 1;
    for i in range(len(contourxs)):
        if contourys[i] < pointy and contourys[j] >= pointy or contourys[j] < pointy and contourys[i] >= pointy:
            if contourxs[i] + (pointy - contourys[i]) / (contourys[j] - contourys[i]) * (contourxs[j] - contourxs[i]) < pointx:
                result = not result
        j = i
    return result

tempTime = time.time()

exportloc = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9fldhnkpreweight"
origfile = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9fldhnkideal\idealdose.dcm"
ptvnames=["PTV_total"]#as for the ptvs, we will try to keep the doses the same between the optimized plan and the ideal.  everywhere else we will minimize.
origrd= pydicom.dcmread(origfile)

for f in listdir(exportloc):
    if isfile(join(exportloc, f)):
        if f[0:2]=="RS":
            rs = pydicom.dcmread(exportloc +"\\"+f)
origdose = origrd.pixel_array * origrd.DoseGridScaling
origin=origrd.ImagePositionPatient
pixelspacing=origrd.PixelSpacing[0]#pixel spacing in the axial plane better be the same x vs y
spacing = origrd.GridFrameOffsetVector[1]-origrd.GridFrameOffsetVector[0]#this is the slice spacing
inside=[]

rois=[]
for ptv in ptvnames:
    for r in rs.StructureSetROISequence:
        if r.ROIName==ptv:
            rois.append(r.ROINumber)

for r in rois:
    for s in rs.ROIContourSequence:
        if s.ReferencedROINumber==r:
            thiscontour = s.ContourSequence
            break
    for k in range(len(thiscontour)):
        c=thiscontour[k].ContourData
        xs = c[::3]
        ys = c[1::3]
        z= np.int32(c[2])#z isnt gonna change per sequence
        minx= np.int32(np.min(xs))
        maxx = np.int32(np.max(xs))
        miny= np.int32(np.min(ys))
        maxy=np.int32(np.max(ys))
        for i in range(minx,maxx):
            for j in range(miny,maxy):
                if IsPointInPolygon(xs,ys,i,j):
                    inside.append([np.round((i-origin[0])/pixelspacing),np.round((j-origin[1])/pixelspacing),np.round((z-origin[2])/spacing)])
inside=np.array(inside, dtype=np.int32)

origdose = origrd.pixel_array * origrd.DoseGridScaling
ptvmask =0*deepcopy(origdose)
for p in inside:
    ptvmask[ p[2] ,p[1],p[0]]=1#whatever this is how the indices are

#first let's get the bounds of the RD array that we want to optimize with.  No need to optimize everywhere
maxorigdose = np.max(origdose)
whererelevant = np.argwhere(origdose>0.5*maxorigdose)#im just lookin at doses within 50% of the max
imin = np.min(whererelevant[:,0])
imax = np.max(whererelevant[:,0])
jmin = np.min(whererelevant[:,1])
jmax = np.max(whererelevant[:,1])
kmin = np.min(whererelevant[:,2])
kmax = np.max(whererelevant[:,2])
goaldose = origdose[imin:imax,jmin:jmax,kmin:kmax]
mask = ptvmask[imin:imax,jmin:jmax,kmin:kmax]#not using this just yet

#get the mlmlc field doses
rds =[]
count=0
for f in listdir(exportloc):
    if isfile(join(exportloc, f)):
        if f[0:2]=="RD":
            tempd = pydicom.dcmread(exportloc +"\\"+f)
            scale = tempd.DoseGridScaling
            temppa = tempd.pixel_array * scale
            rds.append(temppa[imin:imax,jmin:jmax,kmin:kmax])
        if f[0:2]=="RP":
            rp = pydicom.dcmread(exportloc +"\\"+f)
        print("\r loaded file "+ str(count+1) + " of " + str(len(listdir(exportloc))))
        count+=1
#now do the optimization
weights = np.ones(len(rds))
notopt = True
counter = 1
wtd=[0.01,-0.01]#we're gonna search by steps of 1%

allsum = sumall(rds, weights)

#need to first normalize to the mean dose between the mean dose to the sets so that we dont have to search quite as hard.
normalize = np.mean(goaldose)/np.mean(allsum)

print("normalized by " + str(normalize))

weights *= normalize
allsum = sumall(rds, weights)
oldcost = cost(goaldose,allsum)
newcost=oldcost
while notopt:
    notopt=False#if this doesnt switch to True at the deepest part of these loops, exit
    wtdn = deepcopy(wtd)
    for b in range(len(rds)):
        wtdn = deepcopy(wtd)
        otherbeams = deepcopy(rds)
        del otherbeams[b]
        otherweights = deepcopy(weights)
        otherweights= np.delete(otherweights, b)
        othersum = sumall(otherbeams, otherweights) #reduce the number of calculations we have to do when optimizing weight per beam

        while len(wtdn)>0:
            oldcost= newcost
            #weights[b]+=wtdn[0]
            
            if weights[b]+wtdn[0] <= 0:
                oldw=weights[b]
                weights[b]=0
                # if len(wtdn)==1:
                #     del wtdn[0]
            else: weights[b]+=wtdn[0]
            allsum= rds[b]*weights[b] + othersum
            newcost = cost(goaldose,allsum)
            if newcost>oldcost :
                if weights[b]==0: weights[b]=oldw
                    
                else: weights[b] -= wtdn[0]
                
                del wtdn[0]
                newcost=oldcost
            else: 
                print("iteration " + str(counter) + " beam " + str(b+1) + " weight " + str(weights[b])[:6] + " cost " + str(newcost)[:6])
                notopt=True
                if weights[b]==0:
                    del wtdn[0]
    counter+=1
    if counter==100:
        notopt=False#stop after n iterations        
tr.rewriteWeights(origfile, exportloc,weights, "HNMLCfinal")
t = (time.time() - tempTime)
print("Weight optimization completed in " + str(t)[:6]+" seconds")