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

tempTime = time.time()

exportloc = r"K:\Physics Division\Personal folders\ANM\tmp"
origfile = r"C:\Users\anmorrow.BHCS\Desktop\MLC_project\MLMLC_Code\hn178fld\origdose.dcm"
origrd= pydicom.dcmread(origfile)

#first let's get the bounds of the RD array that we want to optimize with.  No need to optimize everywhere
origdose = origrd.pixel_array * origrd.DoseGridScaling
maxorigdose = np.max(origdose)
whererelevant = np.argwhere(origdose>0.5*maxorigdose)#im just lookin at doses within 50% of the max
imin = np.min(whererelevant[:,0])
imax = np.max(whererelevant[:,0])
jmin = np.min(whererelevant[:,1])
jmax = np.max(whererelevant[:,1])
kmin = np.min(whererelevant[:,2])
kmax = np.max(whererelevant[:,2])
goaldose = origdose[imin:imax,jmin:jmax,kmin:kmax]

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
        # if f[0:2]=="RS":
        #     rs = pydicom.dcmread(exportloc +"\\"+f)
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
    if counter==5:
        notopt=False#stop after 3 iterations
#current MU are in rp.FractionGroupSequence[0].ReferencedBeamSequence[b].BeamMeterset.  let's change those out with the new weights.            
tr.rewriteWeights(origfile, exportloc,weights, "HNMLCfinal")
t = (time.time() - tempTime)
print("Weight optimization completed in " + str(t)[:6]+" seconds")