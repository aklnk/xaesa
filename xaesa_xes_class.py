# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:07:33 2018

@author: akali
"""

from numpy import concatenate, logical_and, polyfit, where, zeros, sum
from scipy.integrate import simps
from scipy.interpolate import Rbf

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
        
        self.incidentEnergy = 0
        
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
        
    def removeBackground(self, polPower=1):
        where1 = where( logical_and(self.energy > self.E0, self.energy < self.E1) )
        where2 = where( logical_and(self.energy > self.E2, self.energy < self.E3) )
        
        energy_sliced = concatenate( (self.energy[where1], self.energy[where2]) )
        xes_sliced = concatenate( (self.xes[where1], self.xes[where2]) )

                              
        pf = polyfit(energy_sliced, xes_sliced, polPower)
        
        self.xesBackground = zeros(len(self.energy))
        for i in range(polPower+1):
            self.xesBackground  = self.xesBackground + self.energy**(polPower-i) * pf[i]

        self.xesBkgrCorrected = self.xes - self.xesBackground 
        
    def areaNormalize(self):
        where3 = where( logical_and(self.energy > self.eAreaNormMin, self.energy < self.eAreaNormMax) )    
        
        self.xesAreaNorm = self.xesBkgrCorrected / simps(self.xesBkgrCorrected[where3], self.energy[where3])
        
        print("Integral", simps(self.xesBkgrCorrected[where3], self.energy[where3]))
        print("Sum", sum(self.xesBkgrCorrected[where3]))
        
    def maxNormalize(self):
        pass
    
    def rebin():
        pass      
        
    

class xaesa_transient_class():
    def __init__(self):
        self.name = ""
        self.energy = [] #1D array
        self.transient = []
        self.energySmooth = []
        self.transientSmooth = []
        
        #energy window for smoothing
        self.eSmoothMin = 0
        self.eSmoothMax = 0
        
        self.smoothFactor = 0.0001
        
    def smoothTransient(self, smFactor = 0.0001):
        where3 = where( logical_and(self.energy > self.eSmoothMin, self.energy < self.eSmoothMax ) )
        spl = Rbf( self.energy[where3],self.transient[where3], 
                           function = 'multiquadric', 
                           epsilon = 3, 
                           smooth = smFactor )
        self.energySmooth = self.energy[where3]        
        self.transientSmooth = spl(self.energy[where3])
        
        
class xaesa_kinetic_class():
    def __init__(self):
        self.time = [] #1D array
        self.kinetic = []
        