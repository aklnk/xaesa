# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 11:55:12 2018

@author: akali
"""
import sys
import os
from numpy import   arange, argmin, argmax, asarray, concatenate, copy, delete, \
                    gradient, log, logical_and, multiply, newaxis, pi, power, \
                    sin, cos, sqrt, sum, where, zeros, unique, ones, convolve, isnan

from scipy.optimize import curve_fit
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline, UnivariateSpline
from scipy.signal import symiirorder1

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from .xaesa_constants_formulas import victoreen, windowGauss10, \
                                         me, hbar
                                         
from .xaesa_ft import FT, BFT, GETPHASE, BFTWindow

#import xalglib

class xaesa_exafs_class():
    def __init__(self, init_raw_data_type):
        self.raw_data_type = init_raw_data_type #0 for transmission, 1 for fluorescence, 2 for Mju, 3 for EXAFS
        
        self.name = "" #name
        
        #raw data loaded from experimantal file
        self.energy = [] #1D array
        self.energyRebined = []
        self.energyOriginal = []
        self.i0 = [] #1D array
        self.i1 = [] #1D array
        self.i2 = [] #1D array
        self.ifluo = [] #2D array
        
        #Mju datasets
        self.mju = []
        self.mjuRebined = []
        self.mjuOriginal = []
        self.mjuDerivative= []
        self.victoreen = []
        self.mjuMinusVictoreen = []
        self.mju0 = []

        self.mjuRef = []
        
        #EXAFS datasets
        self.k = []
        self.exafs = []
        self.exafsZeroLine = []
        self.exafsZLC = []
        self.exafsZLCdeglitch = []
        
        # FT datasets
        self.window = []
        self.exafsTimesWindow = []
        self.exafsZLCTimesWindow = []
        self.r = []
        self.fr = []
        self.fi = []        
        self.efr = []
        self.efi = []
        
        self.rZLC = []
        self.frZLC = []
        self.fiZLC = []        
        self.efrZLC = []
        self.efiZLC = []
        
        # BFT datasets
        self.bftWindow = []
        self.bftk = []
        self.bftefr = []
        self.bftefi = []
        self.bftAmp = []
        self.bftPha = []
        self.bftEXAFS = []
        self.bftefrWindow = []
        self.bftefiWindow = []     
        
        #Rebin parameters
        self.dE1 = 3
        self.dE2 = 0.1
        self.dE3 = 1
        
        
        #extraction params
        self.Es = 0
        self.E0 = 0
        self.E1 = 0
        self.E2 = 0
        self.E3 = 0
        
        self.energyShift = 0
        
        self.kPower = 2
        
        self.zeroLineCorr = 0
        
        self.mju0PolinomialDegree = 4
        
        self.normalizationMode = 0 #0 for mju0 normalization, 1 for value normalization at given energy
        self.normalizationEnergy = 0
        
        self.reduceK = 0
        
        self.c = 0
        
        #FT params
        self.kMin = 0.5
        self.kMax = 18
        self.dk = 0.05
        self.rMin = 0
        self.rMax = 6
        self.dr = 0.02
        
        #BFT params
        self.rMinBft = 0
        self.rMaxBft = 6
        self.bftWindowParam = 0.1        
        
        # Deglitching params and result
        self.deglitchIntervals = []
        
        #fit params and results
        self.isFitted = 0
        
        self.fitKMin = 0.5
        self.fitKMax = 15
        self.fitdk = 0.05
        
        self.fitNShels = 1
        
        self.fitParams = []
        self.fitAmps = []
        self.fitPhas = []
        
        self.fitK = []
        self.fitExafs = []
        
        #rdf params and results
        
        self.isRdfed = 0
        
        self.rdfKMin = 0.5
        self.rdfKMax = 15
        self.rdfdk = 0.05
        
        self.rdfRMin = 0.5
        self.rdfRMax = 2.8
        self.rdfdr = 0.01
        
        self.rdfMaxIterations = 20
        
        self.rdfAmpK = []
        self.rdfAmp = []
        self.rdfPhaK = []
        self.rdfPha = []
        
        self.rdfAmpFile = ''
        self.rdfPhaFile = ''
        
        self.rdfK = []
        self.rdfExafs = []
        
        self.rdfR = []
        self.rdf = []
        
    def calculateMjuTransmission(self):
        self.mju = log(self.i0 / self.i1)
        try:
            self.mjuRef = log(self.i1 / self.i2)
        except:
            pass

    
    def calculateMjuFluorescence(self):
        self.mju = self.ifluo.sum(axis=0) / self.i0
        print(self.mju)
        
        
    def estimateEnergyValues(self):
#        try:
#            spl = Rbf( self.energy, self.mju, 
#                       function = 'multiquadric', 
#                       epsilon = 3, #epsilon 3 woks fine
#                       smooth = 5 )
#        except:
#            print("spline error")
#        self.mju = spl(self.energy)
        
        
        #remove 10% of the range if other edges appears at the end of the spectrum
        lastPoint = int(len(self.mju) - len(self.mju)*0.1)
        self.mjuDerivative = gradient(self.mju)
        maxderiv = argmax(self.mjuDerivative[0:lastPoint])
        
        energyMaxDerivValue = self.energy[0:lastPoint][maxderiv]
        
        print("max deriv at ", energyMaxDerivValue)

        self.E0 = self.energy[maxderiv] + 5.5
        self.E1 = self.energy[maxderiv] - (energyMaxDerivValue - self.energy[0])/5
        self.E2 = self.energy[maxderiv] + (self.energy[-1] - energyMaxDerivValue )/5
        self.E3 = self.energy[-1]
        self.Es = self.E1-100
        
        self.kMax = sqrt(  (2*me/hbar**2) * (self.E3-self.E0) * 1.602*10**-19  ) *10**-10 - 0.5
        
    def removeBackground(self):

        if self.Es == 0:
            energyStart = self.energy[0]
        else:
            energyStart = self.Es
        # if energyStart < self.E1-300: #needed if there is another edge before
        #     energyStart = self.E1-300
          
        whereE1 = where( logical_and(self.energy > energyStart, self.energy < self.E1) )
    
        Xb = self.energy[whereE1]
        yb = self.mju[whereE1]
        
        popt, pcov = curve_fit(victoreen, Xb, yb)

        self.victoreen = victoreen(self.energy, popt[0], popt[1])
        self.mjuMinusVictoreen = self.mju - self.victoreen
        
    def findMju0(self):
                
        whereE2E3 = where( logical_and(self.energy > self.E2, self.energy < self.E3) )
    
        energyToFit = self.energy[whereE2E3]
        mjuToFit =  self.mjuMinusVictoreen[whereE2E3]
         
                
        polynomial_features = PolynomialFeatures(self.mju0PolinomialDegree,
                                                 include_bias=True)
        linear_regression = LinearRegression()
    
        pipeline = Pipeline([("polynomial_features", polynomial_features),
                             ("linear_regression", linear_regression)])

        pipeline.fit(energyToFit[:, newaxis], mjuToFit)

        self.mju0 = pipeline.predict(self.energy[:, newaxis])
        
    def calculateEXAFS(self):

        whereE0E3 = where( logical_and(self.energy > self.E0, self.energy < self.E3) )
    
        energyE0E3 = self.energy[whereE0E3]
        mjuE0E3 =  self.mjuMinusVictoreen[whereE0E3]
        mju0E0E3 =  self.mju0[whereE0E3]
        if self.normalizationMode==0: #mju0 normalization
    #        mju0E0E3 =  mju0[whereE0E3]
            mju0E0E3_1 =  self.mju0[whereE0E3]
        if self.normalizationMode==1:
            idx = argmin(abs(self.energy - self.normalizationEnergy))
            mju0E0E3_1 = zeros(len(energyE0E3)) + self.mju0[idx]
        
        self.k = sqrt(  (2*me/hbar**2) * (energyE0E3-self.E0) * 1.602*10**-19  ) *10**-10
    
        self.exafs = (mjuE0E3 - mju0E0E3)/mju0E0E3_1 * power(self.k, self.kPower)
        
        ### test auto reduce points
        if self.reduceK :
            sys.path.append(os.path.join('.', 'lib'))
            import xaslib as xl
            self.k, self.exafs = xl.kPointsReduction(self.k, self.exafs)
        
    def calculateEXAFSZLC(self):
        if self.zeroLineCorr == 0:
            self.exafsZeroLine = zeros( len(self.k))
        if self.zeroLineCorr > 0:  
            ### test code
#            sys.path.append(os.path.join('.', 'lib'))
#            import xaslib as xl
#            nk, ne = xl.kPointsReduction(self.k, self.exafs)
#            spl = Rbf( nk, ne, 
#                       function = 'multiquadric', 
#                       epsilon = 3, #epsilon 3 woks fine
#                       smooth = self.zeroLineCorr )
            ### end test code
            
            spl = Rbf( self.k, self.exafs, 
                       function = 'multiquadric', 
                       epsilon = 3, #epsilon 3 woks fine
                       smooth = self.zeroLineCorr )
            
            
            self.exafsZeroLine = spl(self.k)
        
        self.exafsZLC = self.exafs - self.exafsZeroLine
        
        
    def ftEXAFS(self):

        self.window = windowGauss10(self.k, self.kMin, self.kMax)

        self.exafsTimesWindow = self.exafs * self.window

        self.r, self.fr, self.fi = FT(self.k, self.exafsTimesWindow, self.rMin, self.rMax, self.dr)
        self.efr = zeros(len( self.r))


        sqrt(self.fr*self.fr + self.fi*self.fi, self.efr)
        self.efi = self.fi * (-1)
        
    def ftEXAFSZLC(self):

        self.window = windowGauss10(self.k, self.kMin, self.kMax)

        self.exafsZLCTimesWindow = self.exafsZLC * self.window

        self.rZLC, self.frZLC, self.fiZLC = FT(self.k, self.exafsZLCTimesWindow, self.rMin, self.rMax, self.dr)
        self.efrZLC = zeros(len( self.r))


        sqrt(self.frZLC*self.frZLC + self.fiZLC*self.fiZLC, self.efrZLC)
        self.efiZLC = self.fiZLC * (-1)
        
    def bftEXAFSZLC(self):

        self.bftWindow = BFTWindow(self.rZLC, self.rMinBft, self.rMaxBft, self.bftWindowParam)
        forbftim = self.fiZLC * self.bftWindow
        forbftre = self.frZLC * self.bftWindow
        self.bftefrWindow =  self.efrZLC * self.bftWindow
        self.bftefiWindow =  self.efiZLC * self.bftWindow
        self.bftk,  self.bftefr, self.bftefi =  BFT(self.rZLC, forbftre, forbftim, self.kMin, self.kMax, self.dk)

        wind =  windowGauss10(self.bftk, self.kMin, self.kMax)
        self.bftAmp = sqrt( self.bftefr*self.bftefr / (wind**self.kPower) + 
                                 self.bftefi*self.bftefi / (wind**self.kPower)        )

        self.bftPha = GETPHASE(self.bftefr, self.bftefi)

        self.bftEXAFS = self.bftAmp * sin(self.bftPha)
        
    def changeToRebinedMju(self):
        self.energyOriginal = copy(self.energy)
        self.mjuOriginal = copy(self.mju)
        
        # sort array in there are decreasing points
        sortIndexes = self.energy.argsort()
        sortedEnergy = self.energy[sortIndexes]
        sortedMju = self.mju[sortIndexes]
        uniqueEnergy, uniqueIndexes = unique(sortedEnergy, return_index=True)
    
        spl = InterpolatedUnivariateSpline(uniqueEnergy, sortedMju[uniqueIndexes]) 
        
        kExafsMin = sqrt(  (2*me/hbar**2) * (self.E0-self.E0) * 1.602*10**-19  ) *10**-10
        kExafsMax = sqrt(  (2*me/hbar**2) * (self.E3-self.E0) * 1.602*10**-19  ) *10**-10
        
        kScale = arange(kExafsMin, kExafsMax, 0.025, dtype='float64')
        eScale = kScale**2 * 10**20 / (1.602*10**-19 * (2*me/hbar**2)) + self.E0
        
#        print(kScale, eScale)

        a1 = arange(self.energy[0], self.E1, self.dE1, dtype='float64')
        a2 = arange(self.E1, self.E0, self.dE2, dtype='float64')
#        a3 = arange(self.E2, self.E3, self.dE3)
        self.energyRebined = concatenate((a1, a2, eScale.astype('float64')))
        
        print(self.energyRebined)
        
        self.mjuRebined = spl(self.energyRebined)
        
        self.energy = copy(self.energyRebined)
        self.mju = copy(self.mjuRebined)
        self.mjuDerivative = gradient(self.mju)
        
        self.redoExtraction()
        
    def changeToRebinedMjuAveraging(self, rebinType = 'spline', dE1 = 1, dE2 = 0.1, dK = 0.01, s = 0):
        # rebinType possible values:
        #    'spline'
        #    'rbf-smooth'
        self.energyOriginal = copy(self.energy)
        self.mjuOriginal = copy(self.mju)
        
        # sort array in there are decreasing points
        sortIndexes = self.energy.argsort()
        sortedEnergy = self.energy[sortIndexes]
        sortedMju = self.mju[sortIndexes]
        
        dE = dE1
        newE = []
        newMju = []
        
        #take data before E1
        E1Region = where(  sortedEnergy < self.E1 )
        E1Energy = sortedEnergy[E1Region]
        E1Mju = sortedMju[E1Region]
        newE.append(E1Energy[0])
        newMju.append(E1Mju[0])
        eCenter = E1Energy[0] + dE / 2
        while eCenter + dE/2 < self.E1:
            reg = where( logical_and(sortedEnergy > eCenter-dE/2, sortedEnergy < eCenter+dE/2) )
            valE = sortedEnergy[reg].sum() / len(sortedEnergy[reg])
            val = sortedMju[reg].sum() / len(sortedMju[reg])
            newE.append(valE)
            newMju.append(val)
            eCenter = eCenter + dE
            
        dE = dE2
        
        #take data before E0+50
        E0Region = where(  logical_and(sortedEnergy > self.E1, sortedEnergy < self.E0+50) )
        E0Energy = sortedEnergy[E0Region]
        E0Mju = sortedMju[E0Region]
        newE.append(E0Energy[0])
        newMju.append(E0Mju[0])
        eCenter = E0Energy[0] + dE / 2
        while eCenter + dE/2 < self.E0+50:
            reg = where( logical_and(sortedEnergy > eCenter-dE/2, sortedEnergy < eCenter+dE/2) )
#            print("reg", len(reg), reg)
#            if reg==[]:
#                eCenter = eCenter + dE
            if len(reg[0])==1:
                valE = sortedEnergy[reg][0]
                val = sortedMju[reg][0]
                newE.append(valE)
                newMju.append(val)
            if len(reg[0])>1:   
                valE = sortedEnergy[reg].sum() / len(sortedEnergy[reg])
                val = sortedMju[reg].sum() / len(sortedMju[reg])
                newE.append(valE)
                newMju.append(val)
            eCenter = eCenter + dE
            
            
        dk = dK
        
        #take data after E0
        E3Region = where(  logical_and(sortedEnergy > self.E0+50, sortedEnergy < self.E3) )
        E3Energy = sortedEnergy[E3Region]
        E3Mju = sortedMju[E3Region]
        newE.append(E3Energy[0])
        newMju.append(E3Mju[0])
        
#        kExafsMin = sqrt(  (2*me/hbar**2) * (self.E0+50-self.E0) * 1.602*10**-19  ) *10**-10
#        kExafsMax = sqrt(  (2*me/hbar**2) * (self.E3-self.E0) * 1.602*10**-19  ) *10**-10
        
        kScale = sqrt(  (2*me/hbar**2) * (E3Energy-self.E0) * 1.602*10**-19  ) *10**-10
#        eScale = kScale**2 * 10**20 / (1.602*10**-19 * (2*me/hbar**2)) + self.E0
        kCenter = kScale[0] + dk / 2
        while kCenter + dk/2 < kScale[-1]:
            reg = where( logical_and(kScale > kCenter-dk/2, kScale < kCenter+dk/2) )
#            print("reg", len(reg), reg)
#            if reg==[]:
#                eCenter = eCenter + dE
            if len(reg[0])==1:
                valE = E3Energy[reg][0]
                val = E3Mju[reg][0]
                newE.append(valE)
                newMju.append(val)
            if len(reg[0])>1:   
                valE = E3Energy[reg].sum() / len(E3Energy[reg])
                val = E3Mju[reg].sum() / len(E3Mju[reg])
                newE.append(valE)
                newMju.append(val)
            kCenter = kCenter + dk
            
        
        newE = asarray(newE)
        newMju = asarray(newMju)
        
        kExafsMin = sqrt(  (2*me/hbar**2) * (self.E0+50-self.E0) * 1.602*10**-19  ) *10**-10
        kExafsMax = sqrt(  (2*me/hbar**2) * (self.E3-self.E0) * 1.602*10**-19  ) *10**-10
        
        kScale = arange(kExafsMin, kExafsMax, dk, dtype='float64')
        eScale = kScale**2 * 10**20 / (1.602*10**-19 * (2*me/hbar**2)) + self.E0        

        a1 = arange(self.energy[0], self.E1, self.dE1, dtype='float64')
        a2 = arange(self.E1, self.E0+50, self.dE2, dtype='float64')        
        self.energyRebined = concatenate((a1, a2, eScale.astype('float64')))        
        
        if rebinType == 'spline':
            uniqueEnergy, uniqueIndexes = unique(newE, return_index=True)
            
            newE = newE[uniqueIndexes]
            newMju = newMju[uniqueIndexes]     
            
            newE = newE[~isnan(newE)]
            newMju = newMju[~isnan(newMju)]
        
            spl = InterpolatedUnivariateSpline(newE, newMju, k=1)
            
            
        if rebinType == 'rbf-smooth':
            spl = Rbf( newE, newMju, 
                       function = 'multiquadric', 
                       epsilon = 3, #epsilon 3 woks fine
                       smooth = s )

        
        
        self.mjuRebined = spl(self.energyRebined)
        
        self.mjuRebined = self.mjuRebined[1:-2]
        self.energyRebined = self.energyRebined[1:-2]

        
        self.energy = copy(self.energyRebined)
        self.mju = copy(self.mjuRebined)
        self.mjuDerivative = gradient(self.mju)
        
        
        self.redoExtraction()

        
    def changeToRebinedMjuRbfSmooth(self):
        try:
            spl = Rbf( self.energy, self.mju, 
                       function = 'multiquadric', 
                       epsilon = 3, #epsilon 3 woks fine
                       smooth = 1 )
        except:
            print("rbf error")
        self.mjuRebined = spl(self.energy)
        self.energyRebined = self.energy
        
        self.energy = copy(self.energyRebined)
        self.mju = copy(self.mjuRebined)
        self.mjuDerivative = gradient(self.mju)
      
        self.redoExtraction()
        
    def changeToRebinedMjuAveraging1(self):
        self.energyOriginal = copy(self.energy)
        self.mjuOriginal = copy(self.mju)
        
        # sort array in there are decreasing points
        sortIndexes = self.energy.argsort()
        sortedEnergy = self.energy[sortIndexes]
        sortedMju = self.mju[sortIndexes]
        
        dE = 5
        newE = []
        newMju = []
        
        #take data before E1
        E1Region = where(  sortedEnergy < self.E1 )
        E1Energy = sortedEnergy[E1Region]
        E1Mju = sortedMju[E1Region]
        newE.append(E1Energy[0])
        newMju.append(E1Mju[0])
        eCenter = E1Energy[0] + dE / 2
        while eCenter + dE/2 < self.E1:
            reg = where( logical_and(sortedEnergy > eCenter-dE/2, sortedEnergy < eCenter+dE/2) )
            valE = sortedEnergy[reg].sum() / len(sortedEnergy[reg])
            val = sortedMju[reg].sum() / len(sortedMju[reg])
            newE.append(valE)
            newMju.append(val)
            eCenter = eCenter + dE
            
        dE = 0.1
        
        #take data before E0+50
        E0Region = where(  logical_and(sortedEnergy > self.E1, sortedEnergy < self.E0+50) )
        E0Energy = sortedEnergy[E0Region]
        E0Mju = sortedMju[E0Region]
        newE.append(E0Energy[0])
        newMju.append(E0Mju[0])
        eCenter = E0Energy[0] + dE / 2
        while eCenter + dE/2 < self.E0+50:
            reg = where( logical_and(sortedEnergy > eCenter-dE/2, sortedEnergy < eCenter+dE/2) )
#            print("reg", len(reg), reg)
#            if reg==[]:
#                eCenter = eCenter + dE
            if len(reg[0])==1:
                valE = sortedEnergy[reg][0]
                val = sortedMju[reg][0]
                newE.append(valE)
                newMju.append(val)
            if len(reg[0])>1:   
                valE = sortedEnergy[reg].sum() / len(sortedEnergy[reg])
                val = sortedMju[reg].sum() / len(sortedMju[reg])
                newE.append(valE)
                newMju.append(val)
            eCenter = eCenter + dE
            
            
        dk = 0.01
        
        #take data after E0
        E3Region = where(  logical_and(sortedEnergy > self.E0+50, sortedEnergy < self.E3) )
        E3Energy = sortedEnergy[E3Region]
        E3Mju = sortedMju[E3Region]
        newE.append(E3Energy[0])
        newMju.append(E3Mju[0])
        
#        kExafsMin = sqrt(  (2*me/hbar**2) * (self.E0-self.E0) * 1.602*10**-19  ) *10**-10
#        kExafsMax = sqrt(  (2*me/hbar**2) * (self.E3-self.E0) * 1.602*10**-19  ) *10**-10
        
        kScale = sqrt(  (2*me/hbar**2) * (E3Energy-self.E0) * 1.602*10**-19  ) *10**-10
#        eScale = kScale**2 * 10**20 / (1.602*10**-19 * (2*me/hbar**2)) + self.E0
        kCenter = kScale[0] + dk / 2
        while kCenter + dk/2 < kScale[-1]:
            reg = where( logical_and(kScale > kCenter-dE/2, kScale < kCenter+dE/2) )
#            print("reg", len(reg), reg)
#            if reg==[]:
#                eCenter = eCenter + dE
            if len(reg[0])==1:
                valE = E3Energy[reg][0]
                val = E3Mju[reg][0]
                newE.append(valE)
                newMju.append(val)
            if len(reg[0])>1:   
                valE = E3Energy[reg].sum() / len(E3Energy[reg])
                val = E3Mju[reg].sum() / len(E3Mju[reg])
                newE.append(valE)
                newMju.append(val)
            kCenter = kCenter + dk
            
        
        self.energy = copy(asarray(newE))
        self.mju = copy(asarray(newMju))
        self.energyRebined = copy(asarray(newE))
        self.mjuRebined = copy(asarray(newMju))
        self.mjuDerivative = gradient(self.mju)
        
        self.redoExtraction()
    
    def changeToOriginalMju(self):
        self.energy = copy(self.energyOriginal)
        self.mju = copy(self.mjuOriginal)
        self.mjuDerivative = gradient(self.mju)
        
        self.energyRebined = []
        self.mjuRebined = []
        
        self.redoExtraction()

    def processExpData(self):
        
        if self.raw_data_type == 0: #0 for transmission,
            self.calculateMjuTransmission()
            self.estimateEnergyValues()
            self.removeBackground()
            self.findMju0()
            self.calculateEXAFS()
            self.calculateEXAFSZLC()
            self.ftEXAFS()
            self.ftEXAFSZLC()
            self.bftEXAFSZLC()
            
        if self.raw_data_type == 1: # 1 for fluorescence
            self.calculateMjuFluorescence()
            self.estimateEnergyValues()
            self.removeBackground()
            self.findMju0()
            self.calculateEXAFS()
            self.calculateEXAFSZLC()
            self.ftEXAFS()
            self.ftEXAFSZLC()
            self.bftEXAFSZLC()
            
        if self.raw_data_type == 2: # 2 for Mju
            print("processing", self.raw_data_type)
            self.estimateEnergyValues()
            self.removeBackground()
            self.findMju0()
            self.calculateEXAFS()
            self.calculateEXAFSZLC()
            self.ftEXAFS()
            self.ftEXAFSZLC()
            self.bftEXAFSZLC()
            
        if self.raw_data_type == 3: # 3 for EXAFS
            self.calculateEXAFSZLC()
            self.ftEXAFS()
            self.ftEXAFSZLC()
            self.bftEXAFSZLC()
            
    def changeKPower(self, newKPower):
        actual_k_pow = self.kPower
        new_k_pow = newKPower
        
        self.exafs = self.exafs / self.k ** (actual_k_pow)
        self.exafs = self.exafs * self.k ** (new_k_pow)
        
        self.calculateEXAFSZLC()
        self.ftEXAFS()
        self.ftEXAFSZLC()
        self.bftEXAFSZLC()
        
        self.kPower = newKPower
        
    def redoExtraction(self):
        self.removeBackground()
        self.findMju0()
        self.calculateEXAFS()
        self.calculateEXAFSZLC()
        self.ftEXAFS()
        self.ftEXAFSZLC()
        self.bftEXAFSZLC()
        
    def redoFtBft(self):
        self.ftEXAFS()
        self.ftEXAFSZLC()
        self.bftEXAFSZLC()
        
    def redoBft(self):
        self.bftEXAFSZLC()
        

        
        

            
        



        