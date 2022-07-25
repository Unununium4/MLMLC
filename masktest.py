import pydicom
import numpy as np
from os import listdir
from os.path import isfile, join
import translator as tr
import matplotlib.pyplot as plt
from copy import deepcopy

def IsPointInPolygon(contourxs,contourys, pointx, pointy):
    result = False
    j = len(contourxs) - 1;
    for i in range(len(contourxs)):
        if contourys[i] < pointy and contourys[j] >= pointy or contourys[j] < pointy and contourys[i] >= pointy:
            if contourxs[i] + (pointy - contourys[i]) / (contourys[j] - contourys[i]) * (contourxs[j] - contourxs[i]) < pointx:
                result = not result
        j = i
    return result

exportloc = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9fldhnkpreweight"
origfile = r"K:\Physics Division\Personal folders\ANM\MLMLC\data\9fldhnkideal\idealdose.dcm"

origrd= pydicom.dcmread(origfile)
ptvnames=["PTV_total"]

rds =[]
count=0
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