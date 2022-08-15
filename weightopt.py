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

#calculate time
def timetodeliver(dcmplan,weights, leavesfolder, drate):#drate is in MU/min.  need to include leave position file to get the leaf speeds
    ttd=0#time to deliver
    #assume 600MU/minute.  this is in MU/second - make this an input to allow for FFF rates
    drate = drate/60#convert to MU/sec
    mutimes=[]
    for b in range(len(dcmplan.FractionGroupSequence[0].ReferencedBeamSequence)):
        mutimes.append(dcmplan.FractionGroupSequence[0].ReferencedBeamSequence[b].BeamMeterset * weights[b]/drate)
    gspeed = 6#assume 6 degrees per second gantry rotation (360degree/minute)
    gtimes=[2/gspeed]
    for b in range(1,len(dcmplan.BeamSequence)):
        g1=dcmplan.BeamSequence[b-1].ControlPointSequence[0].GantryAngle
        g2=dcmplan.BeamSequence[b].ControlPointSequence[0].GantryAngle
        gtimes.append(abs(g2-g1)/gspeed)
    lspeed = 25#assume leaf speed of 25mm per second projected at isocenter

    print(len(gtimes))
    print(len(mutimes))

    return ttd

#cost function
def cost(goalptv,goaloar, optptv,optoar): #playing with the cost function.  i think we'll need to just send it points, not the entire 3D sets to speed things up
    maxd = np.max(goalptv)
    n=len(goalptv)+len(goaloar)
    c=0
    for pt in range(len(goalptv)):
        if optptv[pt]< goalptv[pt] or optptv[pt] > maxd:
            c+= abs(optptv[pt]- goalptv[pt])/goalptv[pt]
    for pt in range(len(goaloar)):
        if optoar[pt]>goaloar[pt]:
            c+= (optoar[pt]- goaloar[pt])/goaloar[pt]
    c=100*c/n            
    return c 

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

exportloc = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\H5M5P8_11_22"
origfile = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\HNKRAs\5RADose.dcm"
finalfolder = r"K:\Physics Division\Personal folders\ANM\MLMLC\data"
ptvnames=["PTV_total"]#as for the ptvs, we will try to keep the doses the same between the optimized plan and the ideal.  everywhere else we will minimize.
oarnames=["Parotid_L_P", "Parotid_R_P", "Larynx_P", "Mandible_P","SpinalCanal_P"]
normptv="PTV_6996"
normvol = 95
normrx = 100
rx=69.96
origrd= pydicom.dcmread(origfile)

for f in listdir(exportloc):
    if isfile(join(exportloc, f)):
        if f[0:2]=="RS":
            rs = pydicom.dcmread(exportloc +"\\"+f)
origdose = origrd.pixel_array * origrd.DoseGridScaling
origin=origrd.ImagePositionPatient
pixelspacing=origrd.PixelSpacing[0]#pixel spacing in the axial plane better be the same x vs y
spacing = origrd.GridFrameOffsetVector[1]-origrd.GridFrameOffsetVector[0]#this is the slice spacing
ptvpts=[]
oarpts=[]
normpts = []
#first let's get the bounds of the RD array that we want to optimize with.  No need to optimize everywhere
maxorigdose = np.max(origdose)
whererelevant = np.argwhere(origdose>0.25*maxorigdose)#im just lookin at doses within 25% of the max
imin = np.min(whererelevant[:,0])
imax = np.max(whererelevant[:,0])
jmin = np.min(whererelevant[:,1])
jmax = np.max(whererelevant[:,1])
kmin = np.min(whererelevant[:,2])
kmax = np.max(whererelevant[:,2])

ptvs=[]
oars=[]
for ptv in ptvnames:
    for r in rs.StructureSetROISequence:
        if r.ROIName==ptv:
            ptvs.append(r.ROINumber)
for oar in oarnames:
    for o in rs.StructureSetROISequence:
        if o.ROIName==oar:
            oars.append(o.ROINumber)
for n in rs.StructureSetROISequence:
    if n.ROIName==normptv:
        normnum = n.ROINumber

for r in ptvs:
    for s in rs.ROIContourSequence:
        if s.ReferencedROINumber==r:
            thiscontour = s.ContourSequence
            for k in range(len(thiscontour)):
                c=thiscontour[k].ContourData
                z= c[2]#z isnt gonna change per sequence
                z=np.round((z-origin[2])/spacing)
                if z<imin or z>=imax:
                    continue
                xs = np.multiply(np.add(c[::3],-origin[0]),1/pixelspacing)
                ys = np.multiply(np.add(c[1::3],-origin[1]),1/pixelspacing)
            
                minx= max([np.int32(np.min(xs)), kmin])
                maxx = min([np.int32(np.max(xs)), kmax])
                miny= max([np.int32(np.min(ys)),jmin])
                maxy=min([np.int32(np.max(ys)),jmax])

                for i in range(minx,maxx):
                    for j in range(miny,maxy):
                        if IsPointInPolygon(xs,ys,i,j):
                            #0,1,2-> 2,1,0
                            ptvpts.append([z,j,i])
for o in oars:
    for s in rs.ROIContourSequence:
        if s.ReferencedROINumber==o:
            thiscontour = s.ContourSequence
            for k in range(len(thiscontour)):
                c=thiscontour[k].ContourData
                z= c[2]#z isnt gonna change per sequence
                z=np.round((z-origin[2])/spacing)
                if z<imin or z>=imax:
                    continue
                xs = np.multiply(np.add(c[::3],-origin[0]),1/pixelspacing)
                ys = np.multiply(np.add(c[1::3],-origin[1]),1/pixelspacing)
            
                minx= max([np.int32(np.min(xs)), kmin])
                maxx = min([np.int32(np.max(xs)), kmax])
                miny= max([np.int32(np.min(ys)),jmin])
                maxy=min([np.int32(np.max(ys)),jmax])

                for i in range(minx,maxx):
                    for j in range(miny,maxy):
                        if IsPointInPolygon(xs,ys,i,j):
                            #coordinate axes are 0,1,2-> 2,1,0
                            oarpts.append([z,j,i])

for s in rs.ROIContourSequence:
    if s.ReferencedROINumber==normnum:
        thiscontour = s.ContourSequence
        for k in range(len(thiscontour)):
            c=thiscontour[k].ContourData
            z= c[2]#z isnt gonna change per sequence
            z=np.round((z-origin[2])/spacing)
            if z<imin or z>=imax:
                continue
            xs = np.multiply(np.add(c[::3],-origin[0]),1/pixelspacing)
            ys = np.multiply(np.add(c[1::3],-origin[1]),1/pixelspacing)
        
            minx= max([np.int32(np.min(xs)), kmin])
            maxx = min([np.int32(np.max(xs)), kmax])
            miny= max([np.int32(np.min(ys)),jmin])
            maxy=min([np.int32(np.max(ys)),jmax])

            for i in range(minx,maxx):
                for j in range(miny,maxy):
                    if IsPointInPolygon(xs,ys,i,j):
                        #coordinate axes are 0,1,2-> 2,1,0
                        normpts.append([z,j,i])

oarpts=np.array(oarpts, dtype=np.int32)
ptvpts=np.array(ptvpts, dtype=np.int32)
normpts = np.array(normpts, dtype=np.int32)
imin=min([min(np.transpose(ptvpts)[0]), min(np.transpose(oarpts)[0])])
jmin=min([min(np.transpose(ptvpts)[1]), min(np.transpose(oarpts)[1])])
kmin=min([min(np.transpose(ptvpts)[2]), min(np.transpose(oarpts)[2])])
imax=max([max(np.transpose(ptvpts)[0]), max(np.transpose(oarpts)[0])])+1
jmax=max([max(np.transpose(ptvpts)[1]), max(np.transpose(oarpts)[1])])+1
kmax=max([max(np.transpose(ptvpts)[2]), max(np.transpose(oarpts)[2])])+1
for p in oarpts:
    p[0]+=-imin
    p[1]+=-jmin
    p[2]+=-kmin
for p in ptvpts:
    p[0]+=-imin
    p[1]+=-jmin
    p[2]+=-kmin
for p in normpts:
    p[0]+=-imin
    p[1]+=-jmin
    p[2]+=-kmin

goaldose = origdose[imin:imax,jmin:jmax,kmin:kmax]
goalptvcloud=[]
for pt in ptvpts:
    goalptvcloud.append(goaldose[pt[0],pt[1],pt[2]])
goaloarcloud=[]
for pt in oarpts:
    goaloarcloud.append(goaldose[pt[0],pt[1],pt[2]])

goalptvcloud = np.asarray(goalptvcloud)
goaloarcloud = np.asarray(goaloarcloud)
#get the mlmlc field doses
#rds =[]
preptvcloud=[]
preoarcloud=[]
normcloud = []
count=0
for f in listdir(exportloc):
    if isfile(join(exportloc, f)):
        if f[0:2]=="RD":
            tempd = pydicom.dcmread(exportloc +"\\"+f)
            scale = tempd.DoseGridScaling
            temppa = (tempd.pixel_array * scale)[imin:imax,jmin:jmax,kmin:kmax]
            tempptvcloud=[]
            tempoarcloud=[]
            tempnormcloud=[]
            for pt in ptvpts:
                tempptvcloud.append(temppa[pt[0],pt[1],pt[2]])
            for pt in oarpts:
                tempoarcloud.append(temppa[pt[0],pt[1],pt[2]])
            for pt in normpts:
                tempnormcloud.append(temppa[pt[0],pt[1],pt[2]])
            preptvcloud.append(tempptvcloud)
            preoarcloud.append(tempoarcloud)
            normcloud.append(tempnormcloud)
        if f[0:2]=="RP":
            rp = pydicom.dcmread(exportloc +"\\"+f)
        print("\r loaded file "+ str(count+1) + " of " + str(len(listdir(exportloc))))
        count+=1
#now do the optimization
preptvcloud = np.asarray(preptvcloud)
preoarcloud = np.asarray(preoarcloud)
normcloud = np.asarray(normcloud)
weights = np.ones(len(preptvcloud))
notopt = True
counter = 1
wtd=[0.01,-0.01]#we're gonna search by steps of 1%

#need to first normalize to the mean dose between the mean dose to the sets so that we dont have to search quite as hard.
normalize = np.mean(goalptvcloud)/np.mean(sumall(preptvcloud, weights))

print("normalized by " + str(normalize))

weights *= normalize
allsumptv = sumall(preptvcloud, weights)
alloarptv = sumall(preoarcloud, weights)
oldcost = cost(goalptvcloud,goaloarcloud, allsumptv,alloarptv)
newcost=oldcost
iterchangearray=[]
while notopt:
    notopt=False#if this doesnt switch to True at the deepest part of these loops, exit
    wtdn = deepcopy(wtd)
    itercost= newcost
    for b in range(len(preptvcloud)):
        wtdn = deepcopy(wtd)
        otherbeamsptv = deepcopy(preptvcloud)
        otherbeamsoar = deepcopy(preoarcloud)

        otherbeamsptv=np.delete(otherbeamsptv,b,0)
        otherbeamsoar=np.delete(otherbeamsoar,b,0)

        otherweights = deepcopy(weights)
        otherweights= np.delete(otherweights, b)
        othersumptv = sumall(otherbeamsptv, otherweights) #reduce the number of calculations we have to do when optimizing weight per beam
        othersumoar = sumall(otherbeamsoar, otherweights) 
        while len(wtdn)>0:
            oldcost= newcost
            
            if weights[b]+wtdn[0] <= 0:
                oldw=weights[b]
                weights[b]=0
            else: weights[b]+=wtdn[0]

            allsumptv= preptvcloud[b]*weights[b] + othersumptv
            allsumoar= preoarcloud[b]*weights[b] + othersumoar

            newcost = cost(goalptvcloud,goaloarcloud, allsumptv,allsumoar)
            if newcost>=oldcost :
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
    iterchange = 100*(itercost-newcost)/itercost
    print("Change over this iteration: " + str(iterchange)[:4]+"%")
    iterchangearray.append(iterchange)
    if iterchange <0.1:
        notopt = False#if only 0.1%change in cost from the new iteration go ahead and complete
        print("Optimization completed due to minimal change in cost between iterations.")    
    if counter==100:
        notopt=False#stop after n iterations    
        print("optimization completed due to meeting the iterations limit")   


#renormalize to whatever the original plan was so we get apples to apples
prenorm = sumall(normcloud, weights)
weights = weights* (rx*normrx/100)/np.percentile(prenorm,100-normvol)
#calculate total time and MU for data acquisition


tr.rewriteWeights(origfile, exportloc,weights, finalfolder + r"\H5M5W8_11_22.dcm")



t = (time.time() - tempTime)
np.savetxt(finalfolder+r"\HNMLCcostf.csv",iterchangearray,delimiter =',')



print("Weight optimization completed in " + str(t)[:6]+" seconds")