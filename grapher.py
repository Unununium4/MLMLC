# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:32:05 2019

@author: NickolasT
"""
import matplotlib.pyplot as plt
import numpy as np

# This file converts the fluence data into values that can be printed or graphed
# It also returns summed errors when comparing fluence distributions
# It handles 1D and 2D graphing operations

def mapFluence(x,data1,data2, scale):
    # Graphs the fluence patterns and the difference map
    # It is assumed x is the domain of data1 and data2, (ie their lengths are equal)
    # diffMax defines global difference maximum, make it negative to use local maximum difference
    fig = plt.figure(figsize = (15,5))
    ax = []
    ax.append(fig.add_subplot(1,3,1))
    plt.imshow(data1,vmin=scale[0],vmax=scale[1],interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data1),0])
    plt.title("Optimal Fluence (Interpolated)")
    
    ax.append(fig.add_subplot(1,3,2))
    plt.imshow(data2,vmin=scale[0],vmax=scale[1],interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data2),0])
    plt.title("Approx Fluence through Leaves")
    # plt.colorbar()  #colorbar squishes the 2nd plot
    
    ax.append(fig.add_subplot(1,3,3))
    differenceMap,sumError, sumErrorOver,diffMax = diffMap(data1,data2)
    plt.imshow(differenceMap,interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data1),0])
    plt.title("Relative Difference Map (Blue=Under | Red=Over)")
    return (fig,sumError, sumErrorOver, diffMax)
def mapFluenceREL(x,data1,data2, scale,diffMax,fluenceMax):
    # Graphs the fluence patterns and the difference map
    # It is assumed x is the domain of data1 and data2, (ie their lengths are equal)
    # diffMax defines global difference maximum, make it negative to use local maximum difference
    fig = plt.figure(figsize = (15,5))
    ax = []
    ax.append(fig.add_subplot(1,3,1))
    plt.imshow(np.flip(np.flip(data1,0),1),vmin=scale[0],vmax=scale[1],interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data1),0])
    plt.title("Optimal Fluence (Interpolated)")
    plt.clim(0,fluenceMax)
    
    ax.append(fig.add_subplot(1,3,2))
    plt.imshow(np.flip(np.flip(data2,0),1),vmin=scale[0],vmax=scale[1],interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data2),0])
    plt.title("Approx Fluence through Leaves")
    #plt.colorbar()#colorbar squishes the 2nd plot
    plt.clim(0,fluenceMax)
    
    ax.append(fig.add_subplot(1,3,3))
    differenceMap,sumError, sumErrorOver,diffMax = diffMapREL(data1,data2,diffMax)
    plt.imshow(np.flip(np.flip(differenceMap,0),1),interpolation='nearest', aspect='auto',extent=[min(x),max(x),len(data1),0])
    plt.title("Relative Difference Map (Blue=Under | Red=Over)")
    return (fig,sumError, sumErrorOver, diffMax)
def mapFluence1D(x1,data1,x2,data2):
    # Graphs the fluence pattern for a particular row
    # x1 is the domain of data1 and x2 is the domain of data2
    # Displayed autoscaled graph and the graph scaled to 1
    fig = plt.figure(figsize = (5,5))
    plt.plot(x1,data1)
    plt.plot(x2,data2)
    plt.grid()
    plt.ylim((0,1))
    plt.show()
    return (fig)
    
def convert2RGB(data):
    # *******NOT USED***********
    out = np.zeros([len(data), len(data[0]), 3], dtype=np.uint8)
    tempMax = np.amax(data)
    for i in range(len(data)):
        for j in range(len(data[0])):
            out[i][j][0] = 255*data[i][j]/tempMax
            out[i][j][1] = 255*data[i][j]/tempMax
            out[i][j][2] = 255*data[i][j]/tempMax
    return out

def diffMap(data1, data2):
    # Calculates the difference between corresponding fluence of data1 and data2
    # and normalizes the values. Uses the difference to create RGB values
    # Where white is lowest difference, red is data2 > data1, and blue is data2 < data1
    # RGB scalars are relative to the local maximum of the data set
    
    out = np.zeros([len(data1), len(data1[0]), 3], dtype=np.uint8)
    diffMax = np.amax(np.abs(data2-data1))
    sumError = 0
    sumErrorOver = 0
    for i in range(len(data1)):
        error = 0
        errorOver = 0
        for j in range(len(data1[0])):
            diff = abs(data2[i][j] - data1[i][j])
            error += diff**3
            if (data2[i][j] > data1[i][j]):
                errorOver += diff**3
                out[i][j][0] = 255
                out[i][j][1] = -255*(diff/diffMax)+255
                out[i][j][2] = -255*(diff/diffMax)+255
                
            else:
                out[i][j][0] = -255*(diff/diffMax)+255
                out[i][j][1] = -255*(diff/diffMax)+255
                out[i][j][2] = 255
        sumError += error
        sumErrorOver += errorOver
    
    return out,sumError, sumErrorOver, diffMax
def diffMapREL(data1, data2,diffMax):
    # Calculates the difference between corresponding fluence of data1 and data2
    # and normalizes the values. Uses the difference to create RGB values
    # Where white is lowest difference, red is data2 > data1, and blue is data2 < data1
    # RGB scalars are relative to the maximum specified by diffMax
    
    out = np.zeros([len(data1), len(data1[0]), 3], dtype=np.uint8)
    # If inputted diffMax is negative, then use local diffMax
    sumError = 0
    sumErrorOver = 0
    for i in range(len(data1)):
        error = 0
        errorOver = 0
        for j in range(len(data1[0])):
            diff = abs(data2[i][j] - data1[i][j])
            error += diff**3
            if (data2[i][j] > data1[i][j]):
                errorOver += diff**3
                out[i][j][0] = 255
                out[i][j][1] = -255*(diff/diffMax)+255
                out[i][j][2] = -255*(diff/diffMax)+255
                
            else:
                out[i][j][0] = -255*(diff/diffMax)+255
                out[i][j][1] = -255*(diff/diffMax)+255
                out[i][j][2] = 255
        sumError += error
        sumErrorOver += errorOver
    
    return out,sumError, sumErrorOver, diffMax
def drawLeaves(inSolver,rowEx):
    # Draws a simple diagram of leaf positions, assuming over travel is allowed
    # Draws such that thickest leaf is lowest and thinnest leaf is highest
    # Inputs a solver object and then which row to draw
    yTemp = 0            
    fig = plt.figure()
    for i, mLeaf in enumerate(inSolver.leafLeft[rowEx]):
        X = [inSolver.xMin, mLeaf.pos, mLeaf.pos, inSolver.xMin]
        Y = [yTemp, yTemp, yTemp+(inSolver.nLvls-i),yTemp+(inSolver.nLvls-i)]
        yTemp += (inSolver.nLvls-i)*2
        plt.plot(X,Y)
    yTemp = inSolver.nLvls
    for i, mLeaf in enumerate(inSolver.leafRight[rowEx]):
        X = [inSolver.xMax, mLeaf.pos, mLeaf.pos, inSolver.xMax]
        Y = [yTemp, yTemp, yTemp+(inSolver.nLvls-i),yTemp+(inSolver.nLvls-i)]
        yTemp += (inSolver.nLvls-i)*2-1
        plt.plot(X,Y)
    plt.show()
    inSolver.leafRight[rowEx]
