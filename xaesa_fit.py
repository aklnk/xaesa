# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:55:28 2016

@author: sasha
"""
import os

from .init import QTVer

if QTVer == 4:
    from PyQt4 import QtGui, QtCore
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    
if QTVer == 5:
    from PyQt5 import QtWidgets as QtGui
    from PyQt5 import QtCore
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    
    
import matplotlib.pyplot as plt

import numpy as np

from matplotlib.widgets import SpanSelector

from scipy.interpolate import InterpolatedUnivariateSpline

from scipy.optimize import least_squares, curve_fit

from .xaesa_ft import FT, BFTWindow, BFT, GETPHASE

from .xaesa_constants_formulas import windowGauss10

#from tkinter.filedialog import askopenfilename, asksaveasfilename
#from tkinter import Tk

import gc

def exafsfit_lsq(x, k, exafs, amp, pha, parametrs, var_par, nr_shells, kpow):
    
    #create mixed array with variables and params
    XX = np.zeros(len(var_par))
    varcnt = 0
    parcnt = 0
    for i in range(len(var_par)):
        if(var_par[i]==1): #variable
            XX[i] = x[varcnt]
            varcnt = varcnt + 1
        else:
            XX[i] = parametrs[parcnt]
            parcnt = parcnt + 1
    
#    print("x", x)
    #print("XX", XX)
    chi_model = np.zeros(len(k))
    for i in range(nr_shells):
        chi_model = chi_model + (XX[i*7]/(k*XX[i*7+1]**2)) * amp[i] * \
                    np.exp(-2*XX[i*7+2]*k*k + (2/3)*XX[i*7+4]*k**4 - (4/45)*XX[i*7+6]*k**6) * \
                    np.sin(2*k*XX[i*7+1] - (4/3)*XX[i*7+3]*k**3 + (4/15)*XX[i*7+5]*k**5 + pha[i])
#    chi_model = SO2 * (x[2]/(k*x[3]*x[3])) * amp * exp(-2*x[4]*k*k) * sin(2*k*x[3] + pha)
    return chi_model*k**kpow - exafs

def exafsfit(x, N, R, sigma2):
    
    k = x[0]
    amp = x[1]
    pha= x[2]
    SO2 = x[3]
#    dE0 = X[4]
#    C4 = X[5]
#    C5 = X[6]
#    C6 = X[7]
    chi_model = SO2 * (N/(k*R*R)) * amp * np.exp(-2*sigma2*k*k) * np.sin(2*k*R + pha)
    return chi_model*k*k


class FitWindow(QtGui.QDialog):

    def __init__(self):
        super(FitWindow, self).__init__()
        
        self.bft = []
        self.k = [] 

        self.kamp = [[]]
        self.kpha = [[] ]
        self.amp_orig = [[]]
        self.pha_orig = [[]]

        self.fit_result = []

        self.initUI()

    def initUI(self):
        
        self.shellnr = 1
        self.savedshellnr = 1
        self.isfitted = 0
        
        
        self.fig = plt.figure(3, figsize=(12, 6))
        self.ax_bft = plt.subplot2grid((1,2), (0,0))
        self.ax_bftft = plt.subplot2grid((1,2), (0,1))
        
        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        
        self.fig.tight_layout()
        
#        self.lblNrShells = QtGui.QLabel("Number of shells")
#        self.edtNrShells = QtGui.QLineEdit("1")
        
        self.lblkmin = QtGui.QLabel("K min")
        self.lblkmax = QtGui.QLabel("K max")
        self.lbldk = QtGui.QLabel("dK")
        self.edtkmin = QtGui.QLineEdit("0.5")
        self.edtkmax = QtGui.QLineEdit("15")
        self.edtdk = QtGui.QLineEdit("0.05")
        
        self.lblMaxiterations = QtGui.QLabel("Max number of iterations")
        self.edtMaxiterations = QtGui.QLineEdit("1000")        
        
        self.tabShells = QtGui.QTabWidget()
        self.tabs = []
        self.tabs.append(QtGui.QFrame())
        self.tabShells.addTab(self.tabs[0],"Shell 1")
        
        self.ltShell = []
        self.shellN = []
        self.shellR = []
        self.shellSigma = []
        self.shellC3 = []
        self.shellC4 = []
        self.shellC5 = []
        self.shellC6 = []
#        self.shellE0 = []
        self.shellAmp = []
        self.shellPha = []

        lblN = QtGui.QLabel("N")
        lblR = QtGui.QLabel("R")
        lblSigma = QtGui.QLabel("Sigma")
        lblC3 = QtGui.QLabel("C3")
        lblC4 = QtGui.QLabel("C4")
        lblC5 = QtGui.QLabel("C5")
        lblC6 = QtGui.QLabel("C6")
#        lblE0 = QtGui.QLabel("E0")
        
        lblAmp = QtGui.QLabel("Ampl")
        lblPha = QtGui.QLabel("Phase")

        self.ltShell.append(QtGui.QGridLayout())
        self.shellN.append( [QtGui.QLineEdit("4"), QtGui.QLineEdit("0"), QtGui.QLineEdit("8"), QtGui.QCheckBox()])
        self.shellR.append([QtGui.QLineEdit("2"), QtGui.QLineEdit("0"), QtGui.QLineEdit("4"),  QtGui.QCheckBox()])
        self.shellSigma.append([QtGui.QLineEdit("0.001"), QtGui.QLineEdit("0"), QtGui.QLineEdit("1"),  QtGui.QCheckBox()])
        self.shellC3.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC4.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC5.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC6.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
#        self.shellE0.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellAmp.append(QtGui.QComboBox())
        self.shellPha.append(QtGui.QComboBox())
        
        self.shellAmp[-1].addItem("")
        self.shellPha[-1].addItem("")
        
        self.shellAmp[-1].currentIndexChanged.connect(self.AmpChanged)
        self.shellPha[-1].currentIndexChanged.connect(self.PhaChanged)
    
        
        self.shellN[len(self.shellN)-1][3].setChecked(True)
        self.shellR[len(self.shellR)-1][3].setChecked(True)
        self.shellSigma[len(self.shellSigma)-1][3].setChecked(True)
        
        self.ltShell[0].addWidget(lblN, 0, 0)
        self.ltShell[0].addWidget(lblR, 1, 0)
        self.ltShell[0].addWidget(lblSigma, 2, 0)
        self.ltShell[0].addWidget(lblC3, 3, 0)
        self.ltShell[0].addWidget(lblC4, 4, 0)
        self.ltShell[0].addWidget(lblC5, 5, 0)
        self.ltShell[0].addWidget(lblC6, 6, 0)
#        self.ltShell[0].addWidget(lblE0, 7, 0)  
        self.ltShell[0].addWidget(lblAmp, 7, 0) 
        self.ltShell[0].addWidget(lblPha, 8, 0) 
        
        for i in range(4):           
            self.ltShell[0].addWidget(self.shellN[0][i], 0, 2*i+1)
            self.ltShell[0].addWidget(self.shellR[0][i], 1,  2*i+1)
            self.ltShell[0].addWidget(self.shellSigma[0][i], 2,  2*i+1)
            self.ltShell[0].addWidget(self.shellC3[0][i], 3,  2*i+1)
            self.ltShell[0].addWidget(self.shellC4[0][i], 4,  2*i+1)
            self.ltShell[0].addWidget(self.shellC5[0][i], 5,  2*i+1)
            self.ltShell[0].addWidget(self.shellC6[0][i], 6,  2*i+1)
#            self.ltShell[0].addWidget(self.shellE0[0][i], 7,  2*i+1)
        self.ltShell[0].addWidget(self.shellAmp[0], 7, 1, 1, 7)
        self.ltShell[0].addWidget(self.shellPha[0], 8, 1, 1, 7)
        
        # self.shellAmp[0].addItem("E:/work/development/xaslib/fit/amp0001.dat")
        # self.shellPha[0].addItem("E:/work/development/xaslib/fit/pha0001.dat")
        
        for j in range(7):
            self.ltShell[0].addWidget(QtGui.QLabel("Min. limit"), j, 2)
            self.ltShell[0].addWidget(QtGui.QLabel("Max. limit"), j, 4)
#            self.ltShell[0].addWidget(QtGui.QLabel("Accuracy"), j, 6)
        
        self.tabs[0].setLayout(self.ltShell[0])
        
        self.lblFuncEval = QtGui.QLabel("Number of function evaluations done")
        self.edtFuncEval = QtGui.QLineEdit()
        
        self.lblFitMessage = QtGui.QLabel("Termination reason")
        self.edtFitMessage = QtGui.QLineEdit()
        
        self.lblOptimality = QtGui.QLabel("Cost function")
        self.edtOptimality = QtGui.QLineEdit()
        
        lfit = QtGui.QGridLayout()
        lfit.addWidget(self.lblFitMessage, 0, 0)
        lfit.addWidget(self.edtFitMessage, 0, 1)
        lfit.addWidget(self.lblFuncEval, 1, 0)
        lfit.addWidget(self.edtFuncEval, 1, 1)
        lfit.addWidget(self.lblOptimality, 2, 0)
        lfit.addWidget(self.edtOptimality, 2, 1)
        
        self.btnFit = QtGui.QPushButton('Fit')
        self.btnFit.clicked.connect(self.Fit_leastsquares)

        self.btnSaveFit = QtGui.QPushButton('Save Fit results')
        self.btnSaveFit.clicked.connect(self.saveFit)
        
        self.btnApply = QtGui.QPushButton('Apply')
        self.btnApply.clicked.connect(self.apply)
        
        self.btnCancel = QtGui.QPushButton('Cancel')
        self.btnCancel.clicked.connect(self.cancel)
        
        self.btnAddShell = QtGui.QPushButton('Add shell')
        self.btnAddShell.clicked.connect(self.addshell)
        self.btnRemoveShell = QtGui.QPushButton('Remove shell')
        self.btnRemoveShell.clicked.connect(self.removeshell)
        self.btnOpenAmp = QtGui.QPushButton('Open amplitude file(s) ...')
        self.btnOpenAmp.clicked.connect(self.openamp)
        self.btnOpenPha = QtGui.QPushButton('Open phase file(s) ...')
        self.btnOpenPha.clicked.connect(self.openpha)
        self.btnOpenFeff = QtGui.QPushButton('Open feff file(s) ...')
        self.btnOpenFeff.clicked.connect(self.openfeff)
        
        self.btnSaveFitResults = QtGui.QPushButton('Save fit Results ...')
#        self.btnSaveFitResults.clicked.connect(self.saveFitResults)
        
        lb = QtGui.QGridLayout()
        lb.addWidget(self.btnOpenAmp, 0,0)
        lb.addWidget(self.btnOpenPha, 0,1)
        lb.addWidget(self.btnOpenFeff, 1,0)
        lb.addWidget(self.btnAddShell, 2,0)
        lb.addWidget(self.btnRemoveShell, 2,1)
        
        lfig = QtGui.QGridLayout()
        lfig.addWidget(self.tbar, 0, 0)
        lfig.addWidget(self.canv, 1, 0, 2, 1)
        
        lfig.addLayout(lfit, 3, 0)
        
        lfig.addWidget(self.btnFit, 4, 0)
        lfig.addWidget(self.btnSaveFit, 5, 0)
        lfig.addWidget(self.btnApply, 6, 0)
        lfig.addWidget(self.btnCancel, 7, 0)
        
        lfig.addWidget(self.lblkmin, 0,1)
        lfig.addWidget(self.edtkmin, 0,2)
        lfig.addWidget(self.lblkmax, 0,3)
        lfig.addWidget(self.edtkmax, 0,4)
        lfig.addWidget(self.lbldk, 0,5)
        lfig.addWidget(self.edtdk, 0,6)
        
        lfig.addWidget(self.lblMaxiterations, 1, 1)
        lfig.addWidget(self.edtMaxiterations, 1, 2, 1, 4)
        
        lfig.addWidget(self.tabShells, 2, 1, 2, 6)
        lfig.addLayout(lb, 4,1, 2, 6)
        
        self.setLayout(lfig)
        
    def updateUI(self):
        if self.savedshellnr > 1:
            for i in range(0,self.savedshellnr-1):
                self.addshell()

        self.edtkmin.setText("{:.2f}".format(self.ksettings[0][0]))
        self.edtkmax.setText("{:.2f}".format(self.ksettings[0][1]))
        self.edtdk.setText("{:.2f}".format(self.ksettings[0][2]))
                
        if self.isfitted == 1: #fill with saved fitting params
            self.edtOptimality.setText("{:E}".format(self.costfunction))


            for i in range(self.shellnr):
                self.shellN[i][0].setText("{:.4f}".format(self.fit_params[i][0][0]))
                self.shellN[i][1].setText("{:.4f}".format(self.fit_params[i][0][1]))
                self.shellN[i][2].setText("{:.4f}".format(self.fit_params[i][0][2]))
                self.shellN[i][3].setChecked(bool(self.fit_params[i][0][3]))
                
                self.shellR[i][0].setText("{:.4f}".format(self.fit_params[i][1][0]))
                self.shellR[i][1].setText("{:.4f}".format(self.fit_params[i][1][1]))
                self.shellR[i][2].setText("{:.4f}".format(self.fit_params[i][1][2]))
                self.shellR[i][3].setChecked(bool(self.fit_params[i][1][3]))
                
                self.shellSigma[i][0].setText("{:.4f}".format(self.fit_params[i][2][0]))
                self.shellSigma[i][1].setText("{:.4f}".format(self.fit_params[i][2][1]))
                self.shellSigma[i][2].setText("{:.4f}".format(self.fit_params[i][2][2]))
                self.shellSigma[i][3].setChecked(bool(self.fit_params[i][2][3]))
                
                self.shellC3[i][0].setText("{:.4E}".format(self.fit_params[i][3][0]))
                self.shellC3[i][1].setText("{:.4f}".format(self.fit_params[i][3][1]))
                self.shellC3[i][2].setText("{:.4f}".format(self.fit_params[i][3][2]))
                self.shellC3[i][3].setChecked(bool(self.fit_params[i][3][3]))
                
                self.shellC4[i][0].setText("{:.4E}".format(self.fit_params[i][4][0]))
                self.shellC4[i][1].setText("{:.4f}".format(self.fit_params[i][4][1]))
                self.shellC4[i][2].setText("{:.4f}".format(self.fit_params[i][4][2]))
                self.shellC4[i][3].setChecked(bool(self.fit_params[i][4][3]))
                
                self.shellC5[i][0].setText("{:.4E}".format(self.fit_params[i][5][0]))
                self.shellC5[i][1].setText("{:.4f}".format(self.fit_params[i][5][1]))
                self.shellC5[i][2].setText("{:.4f}".format(self.fit_params[i][5][2]))
                self.shellC5[i][3].setChecked(bool(self.fit_params[i][5][3]))
                
                self.shellC6[i][0].setText("{:.4E}".format(self.fit_params[i][6][0]))
                self.shellC6[i][1].setText("{:.4f}".format(self.fit_params[i][6][1]))
                self.shellC6[i][2].setText("{:.4f}".format(self.fit_params[i][6][2]))
                self.shellC6[i][3].setChecked(bool(self.fit_params[i][6][3]))
                
#            for i in range(int(len(self.fit_amps)/2)):
                self.kamp[i] = self.fit_amps[2*i]
                self.amp_orig[i] = self.fit_amps[2*i+1]
                self.kpha[i] = self.fit_phases[2*i]
                self.pha_orig[i] = self.fit_phases[2*i+1]
#            print(self.fit_amps)

        pass
        
    def Fit_curvefit(self):
        
        kstart = float(self.edtkmin.text())
        kend = float(self.edtkmax.text())
        dk = float(self.edtdk.text())
        common_k = np.arange(kstart, kend, dk)
        
        guess = [0,0,0]
        guess[0] = float(self.shellN[0][0].text())
        guess[1] = float(self.shellR[0][0].text())
        guess[2] = float(self.shellSigma[0][0].text())
        
        varbounds = []
        varbounds.append( ( float(self.shellN[0][1].text()), float(self.shellR[0][1].text()), float(self.shellSigma[0][1].text()) ) )
        varbounds.append( ( float(self.shellN[0][2].text()), float(self.shellR[0][2].text()), float(self.shellSigma[0][2].text()) ) )
        
        kamp, amp_orig = np.genfromtxt("E:/work/development/xaslib/fit/amp0001.dat", usecols=(1,0), unpack=True)
        kpha, pha_orig = np.genfromtxt("E:/work/development/xaslib/fit/pha0001.dat", usecols=(1,0), unpack=True)
        
        splamp = InterpolatedUnivariateSpline(kamp,amp_orig)
        splpha = InterpolatedUnivariateSpline(kpha, pha_orig)
        
        splbft = InterpolatedUnivariateSpline(self.k, self.bft)
        
        amp = splamp(common_k)
        pha = splpha(common_k)
        common_bft = splbft(common_k)
        
#        lsq_result = least_squares(exafsfit, np.array(X), \
#                                   method = 'lm',
##                                   bounds = varbounds,
#                                   args=(self.k, self.bft, amp, pha, 1))
#        print(lsq_result.x)

        x = []
        x.append(common_k)
        x.append(amp)
        x.append(pha)
        x.append(1)

        popt, pcov = curve_fit(exafsfit, x, common_bft , \
                               #method = 'lm',
                               bounds = varbounds,
                               p0 = guess)
        
        self.ax_bft.clear()
        self.ax_bftft.clear()
        self.ax_bft.plot(self.k, self.bft)
#        self.ax_bft.plot(self.k, exafsfit(lsq_result.x, self.k, self.bft, amp, pha, 1)+self.bft)
#        self.ax_bft.plot(self.k, exafsfit(X, self.k, self.bft, amp, pha, 1)+self.bft)
        self.ax_bft.plot(common_k, exafsfit(x, popt[0], popt[1], popt[2]))
        
        print(popt)
        print(pcov)
        
        self.canv.draw()
        
    def Fit_leastsquares(self):
        
        for i in range(self.shellnr):
            if(self.kamp[i]==[]):
                QtGui.QMessageBox.information(self,"Load Amplitude", "Amplitude in shell {:d} not loaded".format(i+1))
                return
            if(self.kpha[i]==[]):
                QtGui.QMessageBox.information(self,"Load Phase", "Phase in shell {:d} not loaded".format(i+1))
                return
        
        kstart = float(self.edtkmin.text())
        kend = float(self.edtkmax.text())
        dk = float(self.edtdk.text())
        self.common_k = np.arange(kstart, kend, dk)
        
        maxiterations = int(self.edtMaxiterations.text())

        #prepare variable and parameter array
        
        splbft = InterpolatedUnivariateSpline(self.k, self.bft)
        self.common_bft = splbft(self.common_k)
        
        varbounds = [[],[]]
        par = []
        var_par = []
        X = []
        edtVarBoxes = []
        amp = []
        pha = []
        for i in range(self.shellnr):
            if self.shellN[i][3].isChecked():
                X.append(float(self.shellN[i][0].text()))
                varbounds[0].append(float(self.shellN[i][1].text()))
                varbounds[1].append(float(self.shellN[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellN[i][0])
            else:
                par.append(float(self.shellN[i][0].text()))
                var_par.append(0)
                
            if self.shellR[i][3].isChecked():
                X.append(float(self.shellR[i][0].text()))
                varbounds[0].append(float(self.shellR[i][1].text()))
                varbounds[1].append(float(self.shellR[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellR[i][0])
            else:
                par.append(float(self.shellR[i][0].text()))
                var_par.append(0)
                
            if self.shellSigma[i][3].isChecked():
                X.append(float(self.shellSigma[i][0].text()))
                varbounds[0].append(float(self.shellSigma[i][1].text()))
                varbounds[1].append(float(self.shellSigma[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellSigma[i][0])
            else:
                par.append(float(self.shellSigma[i][0].text()))
                var_par.append(0)
                
            if self.shellC3[i][3].isChecked():
                X.append(float(self.shellC3[i][0].text()))
                varbounds[0].append(float(self.shellC3[i][1].text()))
                varbounds[1].append(float(self.shellC3[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellC3[i][0])
            else:
                par.append(float(self.shellC3[i][0].text()))
                var_par.append(0)
                
            if self.shellC4[i][3].isChecked():
                X.append(float(self.shellC4[i][0].text()))
                varbounds[0].append(float(self.shellC4[i][1].text()))
                varbounds[1].append(float(self.shellC4[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellC4[i][0])
            else:
                par.append(float(self.shellC4[i][0].text()))
                var_par.append(0)
                
            if self.shellC5[i][3].isChecked():
                X.append(float(self.shellC5[i][0].text()))
                varbounds[0].append(float(self.shellC5[i][1].text()))
                varbounds[1].append(float(self.shellC5[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellC5[i][0])
            else:
                par.append(float(self.shellC5[i][0].text()))
                var_par.append(0)
                
            if self.shellC6[i][3].isChecked():
                X.append(float(self.shellC6[i][0].text()))
                varbounds[0].append(float(self.shellC6[i][1].text()))
                varbounds[1].append(float(self.shellC6[i][2].text()))
                var_par.append(1)
                edtVarBoxes.append(self.shellC6[i][0])
            else:
                par.append(float(self.shellC6[i][0].text()))
                var_par.append(0)
            
            splamp = InterpolatedUnivariateSpline(self.kamp[i], self.amp_orig[i])
            splpha = InterpolatedUnivariateSpline(self.kpha[i], self.pha_orig[i])

            amp.append(splamp(self.common_k))
            pha.append(splpha(self.common_k))
                
            
         
        varbounds[0] = tuple(varbounds[0])
        varbounds[1] = tuple(varbounds[1])
        
        lsq_result = least_squares(exafsfit_lsq, np.array(X), \
#                                   method = 'dogbox',
                                   ftol = 1e-12,  
                                   max_nfev = maxiterations,
                                   bounds = varbounds,
#                                   tr_solver = 'lsmr',
#                                   jac = '3-point',
#                                   loss='soft_l1',
#                                   f_scale=0.1,
                                    verbose = 0,
#                                    x_scale = [1,1,0.001],
                                   args=(self.common_k, self.common_bft, amp, pha, par, var_par, self.shellnr, 2))

        self.edtFuncEval.setText("{:d}".format(lsq_result.nfev))
        self.edtOptimality.setText("{:e}".format(lsq_result.cost))
        self.edtFitMessage.setText(lsq_result.message)

        for i in range(len(lsq_result.x)):
            if i in [0,1,2]:
                edtVarBoxes[i].setText("{:.5f}".format(lsq_result.x[i]))
            else:
                edtVarBoxes[i].setText("{:.2E}".format(lsq_result.x[i]))
        
        self.window = windowGauss10(self.common_k, kstart, kend)
        self.bftw = self.common_bft * self.window
        self.r, self.fr, self.fi = FT(self.common_k, self.bftw, 0, 4, 0.02)
        self.efr = np.sqrt(self.fr*self.fr + self.fi*self.fi)
        self.efi = self.fi * (-1)

        self.fit_result = exafsfit_lsq(lsq_result.x, self.common_k, self.common_bft, amp, pha, par, var_par, self.shellnr, 2)+self.common_bft
        fit_result_w = self.fit_result * self.window
        self.res_r, res_fr, res_fi = FT(self.common_k, fit_result_w, 0, 4, 0.02)
        self.res_efr = np.sqrt(res_fr*res_fr + res_fi*res_fi)
        self.res_efi = res_fi * (-1)
        
        self.ax_bft.clear()
        self.ax_bftft.clear()
        self.ax_bft.plot(self.k, self.bft)
        self.ax_bft.plot(self.common_k, self.fit_result)

        self.ax_bftft.clear()
        
        line1, line2 = self.ax_bftft.plot( self.r,  self.efr, self.r,  self.efi)
        line1.set_color('b')
        line2.set_color('b')
        line2.set_linestyle('dotted')
        
        line1, line2 = self.ax_bftft.plot( self.res_r,  self.res_efr, self.res_r,  self.res_efi)
        line1.set_color('r')
        line2.set_color('r')
        line2.set_linestyle('dotted')
        
        self.canv.draw()

    def apply(self):
        self.fit_params = []
        self.fit_amps = []
        self.fit_phases = []
        self.ksettings = []
        for i in range(self.shellnr):
            self.fit_params.append([ [float(self.shellN[i][0].text()),
                                    float(self.shellN[i][1].text()),
                                    float(self.shellN[i][2].text()),
                                    int(self.shellN[i][3].isChecked())],

                                    [float(self.shellR[i][0].text()),
                                    float(self.shellR[i][1].text()),
                                    float(self.shellR[i][2].text()),
                                    int(self.shellR[i][3].isChecked())],

                                    [float(self.shellSigma[i][0].text()),
                                    float(self.shellSigma[i][1].text()),
                                    float(self.shellSigma[i][2].text()),
                                    int(self.shellSigma[i][3].isChecked())],
            
                                    [float(self.shellC3[i][0].text()),
                                    float(self.shellC3[i][1].text()),
                                    float(self.shellC3[i][2].text()),
                                    int(self.shellC3[i][3].isChecked())],

                                    [float(self.shellC4[i][0].text()),
                                    float(self.shellC4[i][1].text()),
                                    float(self.shellC4[i][2].text()),
                                    int(self.shellC4[i][3].isChecked())],

                                    [float(self.shellC5[i][0].text()),
                                    float(self.shellC5[i][1].text()),
                                    float(self.shellC5[i][2].text()),
                                    int(self.shellC5[i][3].isChecked())],

                                    [float(self.shellC6[i][0].text()),
                                    float(self.shellC6[i][1].text()),
                                    float(self.shellC6[i][2].text()),
                                    int(self.shellC6[i][3].isChecked())]])
            
            self.fit_amps.append( self.kamp[i])
            self.fit_amps.append( self.amp_orig[i])
            self.fit_phases.append( self.kpha[i])
            self.fit_phases.append( self.pha_orig[i])
            
        self.costfunction = float(self.edtOptimality.text())
        
        self.ksettings.append( [float(self.edtkmin.text()),
                                    float(self.edtkmax.text()),
                                    float(self.edtdk.text())] )
                                  
        self.accept()
        
    def cancel(self):
        self.close()
        
    def bftft(self):
        self.window = windowGauss10(self.k, self.k[0], self.k[len(self.k)-1])
        self.bftw = self.bft * self.window
        self.r, self.fr, self.fi = FT(self.k, self.bftw, 0, 4, 0.02)
        self.efr = np.sqrt(self.fr*self.fr + self.fi*self.fi)
        self.efi = self.fi * (-1)
        
        
    def plot(self):
        self.ax_bft.clear()
        self.ax_bftft.clear()
        self.ax_bft.plot(self.k, self.bft)
        
        line1, line2 = self.ax_bftft.plot( self.r,  self.efr, self.r,  self.efi)
        line1.set_color('b')
        line2.set_color('b')
        line2.set_linestyle('dotted')
        
        if(self.fit_result != []):
            kstart = float(self.edtkmin.text())
            kend = float(self.edtkmax.text())
            dk = float(self.edtdk.text())
            self.common_k = np.arange(kstart, kend, dk)
            self.ax_bft.plot(self.common_k, self.fit_result)
        
    def openamp(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)", "Amplitude files (*.amp)"])
#        dlg.setDirectory(os.getcwd())
        if dlg.exec():
            self.fnamp = dlg.selectedFiles()
        else:
            return
        
        self.fnamp.sort()

        for i in range(len(self.shellAmp)):
            self.shellAmp[i].addItems(self.fnamp)
            
    def openpha(self):
        
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)", "Amplitude files (*.pha)"])
#        dlg.setDirectory(os.getcwd())
        if dlg.exec():
            self.fnpha = dlg.selectedFiles()
        else:
            return
        self.fnpha.sort()

        for i in range(len(self.shellPha)):
            self.shellPha[i].addItems(self.fnpha)
            
    def openfeff(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)"])
#        dlg.setDirectory(os.getcwd())
        if dlg.exec():
            self.fnfeff = dlg.selectedFiles()
        else:
            return
        self.fnfeff.sort()
        
        #Extract amplitude and phase from feff files and save to disk
        
        for i in range(len(self.fnfeff)):
            state = 0
            data = []
            f = open(self.fnfeff[i])
            for line in f:
                cols = line.split()
                if cols[0] == '-----------------------------------------------------------------------':
                    state = 1
                    continue
                if cols[0] =='k':
                    state = 2
                    continue
                if state == 1:
                    r = float(cols[2])
                    state = 0
                    continue
                if state == 2:
                    data.append(cols)
        
            new_data_amp = []
            new_data_pha = []
            for j in range(len(data)):
                k = float(data[j][0])
                pha = float(data[j][1]) + float(data[j][3])
                amp = float(data[j][2]) * np.exp( -2 * r / float(data[j][5])) * float(data[j][4])
                new_data_amp.append([k, amp])
                new_data_pha.append([k, pha])

            np.savetxt(self.fnfeff[i] + '.amp', new_data_amp)
            np.savetxt(self.fnfeff[i] + '.pha', new_data_pha)
            
            for j in range(len(self.shellPha)):
                self.shellAmp[j].addItem(self.fnfeff[i] + '.amp')
            for j in range(len(self.shellPha)):
                self.shellPha[j].addItem(self.fnfeff[i] + '.pha')
            
    def addshell(self):
        self.tabs.append(QtGui.QFrame())
        caption = "Shell"+str(self.shellnr+1)
        self.tabShells.addTab(self.tabs[self.shellnr], caption)
        
        lblN = QtGui.QLabel("N")
        lblR = QtGui.QLabel("R")
        lblSigma = QtGui.QLabel("Sigma")
        lblC3 = QtGui.QLabel("C3")
        lblC4 = QtGui.QLabel("C4")
        lblC5 = QtGui.QLabel("C5")
        lblC6 = QtGui.QLabel("C6")
#        lblE0 = QtGui.QLabel("E0")
        
        lblAmp = QtGui.QLabel("Amplitude")
        lblPha = QtGui.QLabel("Phase")

        self.ltShell.append(QtGui.QGridLayout())
        self.shellN.append( [QtGui.QLineEdit("4"), QtGui.QLineEdit("0"), QtGui.QLineEdit("8"),  QtGui.QCheckBox()])
        self.shellR.append([QtGui.QLineEdit("2"), QtGui.QLineEdit("0"), QtGui.QLineEdit("4"),  QtGui.QCheckBox()])
        self.shellSigma.append([QtGui.QLineEdit("0.001"), QtGui.QLineEdit("0"), QtGui.QLineEdit("1"),  QtGui.QCheckBox()])
        self.shellC3.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC4.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC5.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
        self.shellC6.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("-0.1"), QtGui.QLineEdit("0.1"),  QtGui.QCheckBox()])
#        self.shellE0.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellAmp.append(QtGui.QComboBox())
        self.shellPha.append(QtGui.QComboBox())
        
        self.shellAmp[-1].currentIndexChanged.connect(self.AmpChanged)
        self.shellPha[-1].currentIndexChanged.connect(self.PhaChanged)
        
        self.shellN[len(self.shellN)-1][3].setChecked(True)
        self.shellR[len(self.shellR)-1][3].setChecked(True)
        self.shellSigma[len(self.shellSigma)-1][3].setChecked(True)
        
        AllItemsAmp = [self.shellAmp[0].itemText(i) for i in range(self.shellAmp[0].count())]
        AllItemsPha = [self.shellPha[0].itemText(i) for i in range(self.shellPha[0].count())]
        

        self.shellAmp[len(self.shellAmp)-1].addItems(AllItemsAmp)
        self.shellPha[len(self.shellPha)-1].addItems(AllItemsPha)

        
        self.ltShell[self.shellnr].addWidget(lblN, 0, 0)
        self.ltShell[self.shellnr].addWidget(lblR, 1, 0)
        self.ltShell[self.shellnr].addWidget(lblSigma, 2, 0)
        self.ltShell[self.shellnr].addWidget(lblC3, 3, 0)
        self.ltShell[self.shellnr].addWidget(lblC4, 4, 0)
        self.ltShell[self.shellnr].addWidget(lblC5, 5, 0)
        self.ltShell[self.shellnr].addWidget(lblC6, 6, 0)
#        self.ltShell[self.shellnr].addWidget(lblE0, 7, 0)  
        self.ltShell[self.shellnr].addWidget(lblAmp, 7, 0) 
        self.ltShell[self.shellnr].addWidget(lblPha, 8, 0) 
        
        for i in range(4):           
            self.ltShell[self.shellnr].addWidget(self.shellN[self.shellnr][i], 0, 2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellR[self.shellnr][i], 1,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellSigma[self.shellnr][i], 2,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC3[self.shellnr][i], 3,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC4[self.shellnr][i], 4,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC5[self.shellnr][i], 5,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC6[self.shellnr][i], 6,  2*i+1)
#            self.ltShell[self.shellnr].addWidget(self.shellE0[self.shellnr][i], 7,  2*i+1)
        self.ltShell[self.shellnr].addWidget(self.shellAmp[self.shellnr], 7, 1, 1, 7)
        self.ltShell[self.shellnr].addWidget(self.shellPha[self.shellnr], 8, 1, 1, 7)
        

        
        for j in range(7):
            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Min. limit"), j, 2)
            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Max. limit"), j, 4)
#            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Accuracy"), j, 6)
        
        self.tabs[self.shellnr].setLayout(self.ltShell[self.shellnr])
        
        self.kamp.append([])
        self.kpha.append([]) 
        self.amp_orig.append([])
        self.pha_orig.append([])
        
        self.shellnr = self.shellnr +1
        
    def removeshell(self):

        self.tabs.pop()
        
        self.tabShells.removeTab(self.shellnr-1)       


        self.ltShell.pop()
        self.shellN.pop()
        self.shellR.pop()
        self.shellSigma.pop()
        self.shellC3.pop()
        self.shellC4.pop()
        self.shellC5.pop()
        self.shellC6.pop()
#        self.shellE0.pop()
        self.shellAmp.pop()
        self.shellPha.pop()  
        
        self.kamp.pop()
        self.kpha.pop()
        self.amp_orig.pop()
        self.pha_orig.pop()
        
        self.shellnr = self.shellnr -1
        
        gc.collect()
    
    def AmpChanged(self):
        which_shell = -1
        sender = self.sender()
        for i in range(len(self.shellAmp)):
            if self.shellAmp[i] == sender:
                which_shell = i
        if self.shellAmp[which_shell].currentText() == "":
            return      
        ampk, ampo = np.genfromtxt(self.shellAmp[which_shell].currentText(), usecols=(0,1), unpack=True)
        
        self.kamp[which_shell] = ampk
        self.amp_orig[which_shell] = ampo

    def PhaChanged(self):
        which_shell = -1
        sender = self.sender()
        for i in range(len(self.shellPha)):
            if self.shellPha[i] == sender:
                which_shell = i
        if self.shellPha[which_shell].currentText() == "":
            return        
        phak, phao = np.genfromtxt(self.shellPha[which_shell].currentText(), usecols=(0,1), unpack=True)
    
        self.kpha[which_shell] = phak
        self.pha_orig[which_shell] = phao

    def saveFit(self):

        fn = self.savefiledialog_qtgui()
        if fn == "":
            return

        column_captions = ""
        save_data = []
        for i in range(self.shellnr):
            column_captions = column_captions + "Shell{:d} ".format(i)
            values = []
            values.append(float(self.shellN[i][0].text()))
            values.append(float(self.shellR[i][0].text()))
            values.append(float(self.shellSigma[i][0].text()))
            values.append(float(self.shellC3[i][0].text()))
            values.append(float(self.shellC4[i][0].text()))
            values.append(float(self.shellC5[i][0].text()))
            values.append(float(self.shellC6[i][0].text()))

            save_data.append(values)

        np.savetxt(fn + ".fitdata", np.transpose(save_data), header=column_captions)

        column_captions = "k exafs_fit exafs_exp"
        save_array = []
        save_array.append(self.common_k)
        save_array.append(self.fit_result)
        save_array.append(self.common_bft)

        np.savetxt(fn + ".fitexafs", np.transpose(save_array), header=column_captions)

        column_captions = "r_fit ft_real_fit ft_im_fit r_exp ft_real_exp ft_im_exp"
        save_array = []
        save_array.append(self.res_r)
        save_array.append(self.res_efr)
        save_array.append(self.res_efi)
        save_array.append(self.r)
        save_array.append(self.efr)
        save_array.append(self.efi)
        np.savetxt(fn + ".fitft", np.transpose(save_array), header=column_captions)

    def savefiledialog_qtgui(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setAcceptMode(1)  # save dialog
        dlg.setNameFilters(["All files (*.*)"])
        #        dlg.setDirectory(self.currentdir)
        if dlg.exec_():
            flist = dlg.selectedFiles()
            return flist[0]
        else:
            return ""