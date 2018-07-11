# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:07:33 2018

@author: akali
"""

from numpy import concatenate, logical_and, polyfit, where
from scipy.integrate import simps

class xaesa_xes_class():
    def __init__(self):
        
        self.name = ""
        
        #raw data loaded from experimantal file
        self.energy = [] #1D array
        self.energyRebined = []
        self.energyOriginal = []
        
        self.xes = []
        self.xesBackground = []
        self.xesBkgrCorrected = []
        self.xesAreaNorm = []
        self.xesMaxNorm = []
        self.xesRebinned = []
        self.xesOriginal = []
        
        #energy windows for background removal
        self.E0 = 0
        self.E1 = 0
        self.E2 = 0
        self.E3 = 0
        
        #energy window for area normalization
        self.eAreaNormMin = 0
        self.eAreaNormMax = 0
        
        self.eMinRebin = 0
        self.emaxRemin = 0
        self.deRebin = 0
        
    def removeBackground(self):
        where1 = where( logical_and(self.energy > self.E0, self.energy < self.E1) )
        where2 = where( logical_and(self.energy > self.E2, self.energy < self.E3) )
        
        energy_sliced = concatenate( (self.energy[where1], self.energy[where2]) )
        xes_sliced = concatenate( (self.xes[where1], self.xes[where2]) )

                              
        pf = polyfit(energy_sliced, xes_sliced, 1)
        self.xesBackground  = self.energy * pf[0] + pf[1]

        self.xesBkgrCorrected = self.xes - self.xesBackground 
        
    def areaNormalize(self):
        where3 = where( logical_and(self.energy > self.eAreaNormMin, self.energy < self.eAreaNormMax) )    
        
        self.xesAreaNorm = self.xesBkgrCorrected / simps(self.xesBkgrCorrected[where3], self.energy[where3])
        
    def maxNormalize(self):
        pass
    
    def rebin():
        pass
        
        
        
    

class xaesa_transient_class():
    def __init__(self):
        pass