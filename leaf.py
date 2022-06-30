# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:59:02 2019

@author: NickolasT
"""
import numpy as np
from numba.experimental import jitclass
from numba import int32, float32,boolean

spec = [
    ('pos', float32),           # x position of leaf edge
    ('direction',boolean),      # Leaf covers edge (True=Left)
    ('trans', float32),         # Transmission due to leaf attenuation
    ('xMin', float32),          # Minimum x position of leaf
    ('xMax', float32),          # Maximum x position of leaf
    ('level',int32),            # The leaf's level index
    ('row',int32)               # The leaf's row index
]
@jitclass(spec)
class leaf:
    # Object class to store specific information for a leaf
    
    
    def __init__(self,pos, direction, trans, xMin, xMax,level,row):
        self.pos = pos
        self.direction = direction
        self.trans = trans
        self.xMin = xMin
        self.xMax = xMax
        self.level = level
        self.row = row
    
    def mapTrans(self, n):
        # Returns a list of n-elements mapping the resulting tranmission
        # due to the leaf. Where positions covered by the leaf have attenuation
        # and those not covered have no attenuation
        out = np.ones(n)
        step = (self.xMax-self.xMin)/(n-1)
        x = np.arange(self.xMin,self.xMax + step,step)
        if self.direction:
            i = 0
            for i in range(len(x)):
                if (self.pos > x[i]):
                    out[i] = self.trans
            
        else:
            i = n-1
            for i in range(len(x)):
                if(self.pos <= x[i]):
                    out[i] = self.trans
        return out
    
    def calcTrans(self, inFluence):
        # Returns the resulting fluence due to the leaf's attenuation
        return inFluence * self.trans
    def undoTrans(self,inFluence):
        return inFluence / self.trans