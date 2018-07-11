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

from scipy.interpolate import InterpolatedUnivariateSpline

from scipy.optimize import minimize
import scipy
from scipy.integrate import simps

from .xaesa_ft import FT

def exafsrdf_difeval(x, k, exafs, amp_pha_components, dr, kpow):   

    integral = simps(amp_pha_components*x[:,None], dx=dr, axis=0)
    
    sum_sq = sum((integral * k**kpow - exafs)**2)

    return sum_sq

def exafsrdf(x, k, exafs, amp_pha_components, dr, kpow):   
#    integral = simps(transpose(amp_pha_components*x), dx=dr, axis=0)
#    start = timer()
#    integral = simps(amp_pha_components*x[:,None], dx=dr, axis=0)
    integral = simps(amp_pha_components*x[:,None], dx=dr, axis=0)
#    end = timer()
#    print("Time needed:", end - start)
    return integral * k**kpow - exafs

def Window_Gauss10(k, kmin, kmax):
    window = np.zeros(len(k))
    ka = (kmin+kmax)/2.0
    kw = (kmax-kmin)*(kmax-kmin)/9.210340372
    #for i in range(len(k)):
        #window = np.append(window, np.exp(-(k[i]-ka)*(k[i]-ka)/kw))
    knp = np.asarray(k, float)
    wp = -(knp-ka)*(knp-ka)/kw
    np.exp(wp, window)
        
    return window

def Gaussian(rdf_r, A, xc, w):
    return A * np.exp( -(rdf_r-xc)**2 / (2*w) )
    
class RDFinder(QtCore.QObject):
    rdf = []
    options1 = ()
    varbounds1 = []
    args = ()
    
    updateProgress  = QtCore.pyqtSignal(np.ndarray) 
    finished = QtCore.pyqtSignal(scipy.optimize.OptimizeResult)
    
    def __init__(self, parent=None):
        QtCore.QObject.__init__(self, parent)
        print("Init RDFinder")
        
    def calculate(self):
        print("start calculation")
        rdf_result = minimize(exafsrdf_difeval, self.rdf, 

                                    method='L-BFGS-B',
#                                    method='TNC', #BAD
#                                    method='SLSQP',
                                    bounds = self.varbounds1,

                                    callback = self.callbackF,
                                    options = self.options1,
                                    args=self.args ) 
        
        self.finished.emit(rdf_result)
        
    def callbackF(self, Xi):
        self.updateProgress.emit(Xi)
        
        pass
    
class MyStream(QtCore.QObject):
    message = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(MyStream, self).__init__(parent)

    def write(self, message):
        self.message.emit(str(message))

class RdfWindow(QtGui.QDialog):

    def __init__(self):
        super(RdfWindow, self).__init__()
        
        self.bft = []
        self.k = []   
        self.kpow = 0  
#        self.rdf_r = []
#        self.rdf = []
        self.params = []
        self.amp_orig = [[],[]] #0 - k, 1 - amp or pha
        self.pha_orig = [[],[]]
        self.rdf_exafs = []  
        self.rdf_exafs_k = []        

        self.initUI()

    @QtCore.pyqtSlot(str)
    def on_myStream_message(self, message):
        self.textEdit.moveCursor(QtGui.QTextCursor.End)
        self.textEdit.insertPlainText(message)

    def initUI(self):        
        
        self.fig = plt.figure(5, figsize=(18, 6))
        self.ax_bft = plt.subplot2grid((1,3), (0,0))
        self.ax_bftft = plt.subplot2grid((1,3), (0,1))
        self.ax_rdf = plt.subplot2grid((1,3), (0,2))
        
        
        self.ax_bft.set_xlabel('Wavevector k, $\AA^{-1}$')
        self.ax_bft.set_ylabel('EXAFS, $\AA^{-2}$')
        self.ax_bftft.set_xlabel('Distance R, $\AA$')
        self.ax_bftft.set_ylabel('Fourier transform, $\AA^{-3}$')
        self.ax_rdf.set_xlabel('Distance R, $\AA$')
        self.ax_rdf.set_ylabel('RDF G(R)')
        
        self.line1_bft, = self.ax_bft.plot( [],  [], label = 'BFT' )
        self.line2_bft, = self.ax_bft.plot( [],  [], label = 'fit' )
        self.line1_bft.set_color('b')
        self.line2_bft.set_color('g')
        self.ax_bft.legend(handles=[self.line1_bft, self.line2_bft])
        
        self.line1_bftft, self.line2_bftft = self.ax_bftft.plot( [],  [], [],  [])
        self.line1_bftft.set_color('b')
        self.line2_bftft.set_color('b')
        self.line2_bftft.set_linestyle('dotted')
        
        self.line3_bftft, self.line4_bftft = self.ax_bftft.plot( [],  [], [],  [])
        self.line3_bftft.set_color('g')
        self.line4_bftft.set_color('g')
        self.line4_bftft.set_linestyle('dotted')
        
        self.line1_rdf, = self.ax_rdf.plot( [],  [], marker='o' )
        self.line2_rdf, = self.ax_rdf.plot( [],  [] )
        
        self.fig.tight_layout()
        
        
        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        

        self.rdf_r = np.arange(0.8, 2.5, 0.01)
        
        self.rdf = np.zeros(len(self.rdf_r))

        
#        self.lblNrShells = QtGui.QLabel("Number of shells")
#        self.edtNrShells = QtGui.QLineEdit("1")
        
        self.lblkmin = QtGui.QLabel("K min")
        self.lblkmax = QtGui.QLabel("K max")
        self.lbldk = QtGui.QLabel("dK")
        self.lblrmin = QtGui.QLabel("Rmin")
        self.lblrmax = QtGui.QLabel("Rmax")
        self.lbldr = QtGui.QLabel("dR")
        self.edtkmin = QtGui.QLineEdit("0.5")
        self.edtkmax = QtGui.QLineEdit("15")
        self.edtdk = QtGui.QLineEdit("0.05")
        self.edtrmin = QtGui.QLineEdit("0.5")
        self.edtrmax = QtGui.QLineEdit("2.8")
        self.edtdr = QtGui.QLineEdit("0.01")
        
        self.edtrmin.editingFinished.connect(self.resetFit)
        self.edtrmax.editingFinished.connect(self.resetFit)
        self.edtdr.editingFinished.connect(self.resetFit)
        
#        self.lblCentralDist = QtGui.QLabel("Central distance")
#        self.edtCentralDist = QtGui.QLineEdit("1.8")
#        self.lblWidth = QtGui.QLabel("Width")
#        self.edtWidth = QtGui.QLineEdit("0.001")
#        self.lblAmplitude = QtGui.QLabel("Amplitude")
#        self.edtAmplitude = QtGui.QLineEdit("50")
        
        self.lblIterations = QtGui.QLabel("Max. number of iterations")
        self.edtIterations = QtGui.QLineEdit("20")
        self.chkFirstFit = QtGui.QCheckBox("Start new fit (checked) / continue fit (unchecked)")
        self.chkFirstFit.setChecked(True)
        
        self.Amp = QtGui.QComboBox()
        self.Pha = QtGui.QComboBox()
        
        lblAmp = QtGui.QLabel("Amplitude")
        lblPha = QtGui.QLabel("Phase")
        
        self.Amp.addItem("")
        self.Amp.addItem("E:/work/development/xaslib/fit/amp0001.dat")
        self.Pha.addItem("")
        self.Pha.addItem("E:/work/development/xaslib/fit/pha0001.dat")
        
        self.lblFuncEval = QtGui.QLabel("Number of function evaluations done")
        self.edtFuncEval = QtGui.QLineEdit()
        
        self.lblFitMessage = QtGui.QLabel("Iterations done")
        self.edtFitMessage = QtGui.QLineEdit()
        
        self.lblOptimality = QtGui.QLabel("Residual")
        self.edtOptimality = QtGui.QLineEdit()
        
        lfit = QtGui.QGridLayout()
        lfit.addWidget(self.lblFitMessage, 0, 0)
        lfit.addWidget(self.edtFitMessage, 0, 1)
        lfit.addWidget(self.lblFuncEval, 1, 0)
        lfit.addWidget(self.edtFuncEval, 1, 1)
        lfit.addWidget(self.lblOptimality, 2, 0)
        lfit.addWidget(self.edtOptimality, 2, 1)
        
        self.btnFit = QtGui.QPushButton('Fit')
        self.btnFit.clicked.connect(self.Fit_rdf)
        
        self.btnSaveRdf = QtGui.QPushButton('Save RDF')
        self.btnSaveRdf.clicked.connect(self.saveRDF)
        
        self.btnApply = QtGui.QPushButton('Apply')
        self.btnApply.clicked.connect(self.apply)
        
        self.btnCancel = QtGui.QPushButton('Cancel')
        self.btnCancel.clicked.connect(self.cancel)
        
        self.btnOpenAmp = QtGui.QPushButton('Open amplitude file(s) ...')
        self.btnOpenAmp.clicked.connect(self.openamp)
        self.btnOpenPha = QtGui.QPushButton('Open phase file(s) ...')
        self.btnOpenPha.clicked.connect(self.openpha)
        self.btnOpenFeff = QtGui.QPushButton('Open feff file(s) ...')
        self.btnOpenFeff.clicked.connect(self.openfeff)
        
#        self.textEdit = QtGui.QTextEdit(self)
        
        
#        lb.addWidget(self.textEdit, 1, 0)
        
        lfig = QtGui.QGridLayout()
        lfig.addWidget(self.tbar, 0, 0, 1, 2)
        lfig.addWidget(self.canv, 1, 0, 1, 2)
        
        lfig.addLayout(lfit, 2, 0)
        
        lfig.addWidget(self.btnFit, 3, 0)
        lfig.addWidget(self.btnSaveRdf, 4, 0)
        lfig.addWidget(self.btnApply, 5, 0)
        lfig.addWidget(self.btnCancel, 6, 0)
        
        lp = QtGui.QGridLayout()
        lp.addWidget(self.lblkmin, 0,0)
        lp.addWidget(self.edtkmin, 0,1)
        lp.addWidget(self.lblkmax, 0,2)
        lp.addWidget(self.edtkmax, 0,3)
        lp.addWidget(self.lbldk, 0,4)
        lp.addWidget(self.edtdk, 0,5)
        
        lp.addWidget(self.lblrmin, 1,0)
        lp.addWidget(self.edtrmin, 1,1)
        lp.addWidget(self.lblrmax, 1,2)
        lp.addWidget(self.edtrmax, 1,3)
        lp.addWidget(self.lbldr, 1,4)
        lp.addWidget(self.edtdr, 1,5)
        
#        lfig.addWidget(self.lblAmplitude, 2,1)
#        lfig.addWidget(self.edtAmplitude, 2,2)
#        lfig.addWidget(self.lblCentralDist, 2,3)
#        lfig.addWidget(self.edtCentralDist, 2,4)
#        lfig.addWidget(self.lblWidth, 2,5)
#        lfig.addWidget(self.edtWidth, 2,6)
        
        lp.addWidget(self.lblIterations, 2, 0, 1, 3)
        lp.addWidget(self.edtIterations, 2, 3, 1, 3)
        lp.addWidget(self.chkFirstFit, 3, 0, 1, 6)
        
        lp.addWidget(lblAmp, 4,0)
        lp.addWidget(self.Amp, 4,1, 1, 4)
        lp.addWidget(lblPha, 5,0)
        lp.addWidget(self.Pha, 5,1, 1, 4)
        lp.addWidget(self.btnOpenAmp, 4, 5)
        lp.addWidget(self.btnOpenPha, 5, 5)
        lp.addWidget(self.btnOpenFeff, 6, 0)
        
        lfig.addLayout(lp, 2,1, 5, 1)

        
        
        self.setLayout(lfig)
        
        self.canv.draw()
        
        
    def Fit_rdf(self):
        
        self.kstart = float(self.edtkmin.text())
        self.kend = float(self.edtkmax.text())
        self.dk = float(self.edtdk.text())
        self.common_k = np.arange(self.kstart, self.kend, self.dk)
        
        self.rstart = float(self.edtrmin.text())
        self.rend = float(self.edtrmax.text())
        self.dr = float(self.edtdr.text())         
        
        self.rdf_r = np.arange(self.rstart, self.rend, self.dr)

        if self.chkFirstFit.isChecked():
            self.rdf = np.zeros(len(self.rdf_r))
    

        #prepare variable and parameter array
        
        splbft = InterpolatedUnivariateSpline(self.k, self.bft)
        self.common_bft = splbft(self.common_k)
        
        varbounds = [[],[]]
        varbounds1 = []

        for i in range(len(self.rdf_r)):
            varbounds[0].append(0)
            varbounds[1].append(1E6)
            varbounds1.append((0,1E6))
        varbounds[0] = tuple(varbounds[0])
        varbounds[1] = tuple(varbounds[1])

        par = []
        var_par = []
        edtVarBoxes = []     

        if  self.amp_orig[0] == [] or self.pha_orig[0] == []:
                
            kamp, amp_orig = np.genfromtxt(self.Amp.currentText(), usecols=(0,1), unpack=True)
            kpha, pha_orig = np.genfromtxt(self.Pha.currentText(), usecols=(0,1), unpack=True)
            self.amp_orig[0] = kamp
            self.pha_orig[0] = kpha
            self.amp_orig[1] = amp_orig
            self.pha_orig[1] = pha_orig
        else:
            kamp = self.amp_orig[0]
            amp_orig = self.amp_orig[1]
            kpha = self.pha_orig[0]
            pha_orig = self.pha_orig[1]
        
        print(kamp, amp_orig)
        splamp = InterpolatedUnivariateSpline(kamp, amp_orig)
        splpha = InterpolatedUnivariateSpline(kpha, pha_orig)

        amp = splamp(self.common_k)
        pha = splpha(self.common_k)            
        
        pha_component = np.sin(2*np.transpose(np.array([self.rdf_r])).dot(np.array([self.common_k])) + pha)

        rdf_r2 = self.rdf_r*self.rdf_r

        amp_component = amp * (1/ np.transpose(np.array([rdf_r2])).dot(np.array([self.common_k])) )

        self.amptimespha = amp_component * pha_component

        options1={'disp': False, 
                 'maxls': 20, 
                 'iprint': -1, 
                 'gtol': 1e-05, 
                 'eps': 1e-08, 
                 'maxiter': int(self.edtIterations.text()), 
                 'ftol': 2.220446049250313e-09, 
                 'maxcor': 10, 
                 'maxfun': 15000}
        
        args=(self.common_k, self.common_bft, self.amptimespha, self.dr, self.kpow)

        self.rdffnd = RDFinder()
        self.rdffnd.args = args
        self.rdffnd.options1 = options1
        self.rdffnd.varbounds1 = varbounds1
        self.rdffnd.rdf = self.rdf
        
        self.wthread = QtCore.QThread(self)        
        self.rdffnd.moveToThread(self.wthread)
        self.wthread.started.connect(self.rdffnd.calculate)
        self.rdffnd.finished.connect(self.RDFinished)
        self.rdffnd.updateProgress.connect(self.updateProgress)
        self.rdffnd.finished.connect(self.wthread.quit)
        print("Starting optimisation...")
        self.wthread.start()
        print("optimisation started...")


    def updateProgress(self, Xi):
        fresidual = exafsrdf_difeval(Xi, self.common_k, self.common_bft, self.amptimespha, self.dr, self.kpow)
        print("Step residual: ", fresidual)
        self.edtOptimality.setText("{:e}".format(fresidual))
        
        self.window = Window_Gauss10(self.common_k, self.kstart, self.kend)
        self.bftw = self.common_bft * self.window
        self.r, self.fr, self.fi = FT(self.common_k, self.bftw, 0, 4, 0.02)
        self.efr = np.sqrt(self.fr*self.fr + self.fi*self.fi)
        self.efi = self.fi * (-1)
        
        fit_result = exafsrdf(Xi, self.common_k, self.common_bft, self.amptimespha, self.dr, self.kpow)+self.common_bft
        fit_result_w = fit_result * self.window
        res_r, res_fr, res_fi = FT(self.common_k, fit_result_w, 0, 4, 0.02)
        res_efr = np.sqrt(res_fr*res_fr + res_fi*res_fi)
        res_efi = res_fi * (-1)
        
        self.line1_bft.set_xdata(self.k)
        self.line1_bft.set_ydata(self.bft)
        self.line2_bft.set_xdata(self.common_k)
        self.line2_bft.set_ydata(fit_result)
        
        self.line1_bftft.set_xdata(self.r)
        self.line1_bftft.set_ydata(self.efr)
        
        self.line2_bftft.set_xdata(self.r)
        self.line2_bftft.set_ydata(self.efi)
        
        self.line3_bftft.set_xdata(res_r)
        self.line3_bftft.set_ydata(res_efr)
        
        self.line4_bftft.set_xdata(res_r)
        self.line4_bftft.set_ydata(res_efi)
        
        self.line1_rdf.set_xdata(self.rdf_r)
        self.line1_rdf.set_ydata(Xi)
        self.ax_bft.relim()
        self.ax_bft.autoscale()
        self.ax_bftft.relim()
        self.ax_bftft.autoscale()
        self.ax_rdf.relim()
        self.ax_rdf.autoscale()
        
        self.canv.draw()
        
    def RDFinished(self, lsq_result):
        
        print("RDF finnished")
        
        optim = exafsrdf_difeval(lsq_result.x, self.common_k, self.common_bft, self.amptimespha, self.dr, self.kpow)

        self.edtFuncEval.setText("{:d}".format(lsq_result.nfev))
        self.edtOptimality.setText("{:e}".format(optim))
        self.edtFitMessage.setText("{:d}".format(lsq_result.nit))        
        
        self.window = Window_Gauss10(self.common_k, self.kstart, self.kend)
        self.bftw = self.common_bft * self.window
        self.r, self.fr, self.fi = FT(self.common_k, self.bftw, 0, 4, 0.02)
        self.efr = np.sqrt(self.fr*self.fr + self.fi*self.fi)
        self.efi = self.fi * (-1)
        
        self.fit_result = exafsrdf(lsq_result.x, self.common_k, self.common_bft, self.amptimespha, self.dr, self.kpow)+self.common_bft
        fit_result_w = self.fit_result * self.window
        res_r, res_fr, res_fi = FT(self.common_k, fit_result_w, 0, 4, 0.02)
        res_efr = np.sqrt(res_fr*res_fr + res_fi*res_fi)
        res_efi = res_fi * (-1)
        
        self.rdf_exafs = self.fit_result
        self.rdf_exafs_k = self.common_k
        
        self.line1_bft.set_xdata(self.k)
        self.line1_bft.set_ydata(self.bft)
        self.line2_bft.set_xdata(self.common_k)
        self.line2_bft.set_ydata(self.fit_result)
        
        self.line1_bftft.set_xdata(self.r)
        self.line1_bftft.set_ydata(self.efr)
        
        self.line2_bftft.set_xdata(self.r)
        self.line2_bftft.set_ydata(self.efi)
        
        self.line3_bftft.set_xdata(res_r)
        self.line3_bftft.set_ydata(res_efr)
        
        self.line4_bftft.set_xdata(res_r)
        self.line4_bftft.set_ydata(res_efi)
        
        self.line1_rdf.set_xdata(self.rdf_r)
        self.line1_rdf.set_ydata(lsq_result.x)
        
        self.line2_rdf.set_xdata(self.rdf_r)
        self.line2_rdf.set_ydata(self.rdf)
        
        
        self.ax_bft.relim()
        self.ax_bft.autoscale()
        self.ax_bftft.relim()
        self.ax_bftft.autoscale()
        self.ax_rdf.relim()
        self.ax_rdf.autoscale()
        
        self.canv.draw()
        
        self.rdf = lsq_result.x
        self.chkFirstFit.setChecked(False)
        
    def saveRDF(self):
        fn = self.savefiledialog_qtgui()
        if fn == "":
            return 
        column_captions = "#r rdf k"
        save_array = []
        save_array.append(self.rdf_r)
        save_array.append(self.rdf)        
        np.savetxt(fn+".rdf", np.transpose(save_array), header=column_captions)
        
        column_captions = "#k exafs_rdf"
        save_array = []
        save_array.append(self.common_k)
        save_array.append(self.fit_result)
        np.savetxt(fn+".rdfexafs", np.transpose(save_array), header=column_captions)
        
        column_captions = "#r ft_real ft_im"
        save_array = []
        save_array.append(self.r)
        save_array.append(self.efr)
        save_array.append(self.efi)
        np.savetxt(fn+".rdfft", np.transpose(save_array), header=column_captions)
        
    def apply(self):
        self.params = []
        self.params.append(float(self.edtkmin.text()))
        self.params.append(float(self.edtkmax.text()))
        self.params.append(float(self.edtdk.text()))
        self.params.append(float(self.edtrmin.text()))
        self.params.append(float(self.edtrmax.text()))
        self.params.append(float(self.edtdr.text()))
        self.params.append(float(self.edtIterations.text()))
        self.accept()
        
    def cancel(self):
        self.close()
        
    def bftft(self):
        self.window = Window_Gauss10(self.k, self.k[0], self.k[len(self.k)-1])
        self.bftw = self.bft * self.window
        self.r, self.fr, self.fi = FT(self.k, self.bftw, 0, 4, 0.02)
        self.efr = np.sqrt(self.fr*self.fr + self.fi*self.fi)
        self.efi = self.fi * (-1)
        
        
    def plot(self):
        
        self.line1_bft.set_xdata(self.k)
        self.line1_bft.set_ydata(self.bft)
        
        self.line2_bft.set_xdata(self.rdf_exafs_k)
        self.line2_bft.set_ydata(self.rdf_exafs)
        
        self.line1_bftft.set_xdata(self.r)
        self.line1_bftft.set_ydata(self.efr)
        
        self.line2_bftft.set_xdata(self.r)
        self.line2_bftft.set_ydata(self.efi)
        
        self.line1_rdf.set_xdata(self.rdf_r)
        self.line1_rdf.set_ydata(self.rdf)
        
        print(len(self.rdf_exafs_k))
        print(len(self.rdf_exafs))
        
        self.ax_bft.relim()
        self.ax_bft.autoscale()
        self.ax_bftft.relim()
        self.ax_bftft.autoscale()
        self.ax_rdf.relim()
        self.ax_rdf.autoscale()
        
        self.canv.draw()
        
    def openamp(self):

        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)", "Amplitude files (*.amp)"])
        dlg.setDirectory(os.getcwd())
        if dlg.exec_():
            self.fnamp = dlg.selectedFiles()
        else:
            return

        self.fnamp.sort()

        self.Amp.clear()
        self.Amp.addItems(self.fnamp)
            
    def openpha(self):

        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)", "Amplitude files (*.pha)"])
        dlg.setDirectory(os.getcwd())
        if dlg.exec_():
            self.fnpha = dlg.selectedFiles()
        else:
            return
        self.fnpha.sort()

        self.Pha.clear()
        self.Pha.addItems(self.fnpha)
        
    def openfeff(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["All files (*.*)"])
        dlg.setDirectory(os.getcwd())
        if dlg.exec_():
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
            
#            print(new_data_amp)
#            print(new_data_pha)
            np.savetxt(self.fnfeff[i] + '.amp', new_data_amp)
            np.savetxt(self.fnfeff[i] + '.pha', new_data_pha)
            
            self.Amp.clear()
            self.Amp.addItem(self.fnfeff[i] + '.amp')
            
            self.Pha.clear()
            self.Pha.addItem(self.fnfeff[i] + '.pha')
            
            
            
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
        lblE0 = QtGui.QLabel("E0")
        
        lblAmp = QtGui.QLabel("Amplitude")
        lblPha = QtGui.QLabel("Phase")

        self.ltShell.append(QtGui.QGridLayout())
        self.shellN.append( [QtGui.QLineEdit("4"), QtGui.QLineEdit("0"), QtGui.QLineEdit("8"), QtGui.QLineEdit("0.00001"), QtGui.QCheckBox()])
        self.shellR.append([QtGui.QLineEdit("2"), QtGui.QLineEdit("0"), QtGui.QLineEdit("4"), QtGui.QLineEdit("0.0000001"), QtGui.QCheckBox()])
        self.shellSigma.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("1"), QtGui.QLineEdit("0.00000001"), QtGui.QCheckBox()])
        self.shellC3.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellC4.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellC5.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellC6.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellE0.append([QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0"), QtGui.QLineEdit("0.0001"), QtGui.QCheckBox()])
        self.shellAmp.append(QtGui.QComboBox())
        self.shellPha.append(QtGui.QComboBox())
        
        self.shellAmp[len(self.shellAmp)-1].addItem("E:/work/development/xaslib/fit/amp0001.dat")
        self.shellPha[len(self.shellPha)-1].addItem("E:/work/development/xaslib/fit/pha0001.dat")
        
        self.ltShell[self.shellnr].addWidget(lblN, 0, 0)
        self.ltShell[self.shellnr].addWidget(lblR, 1, 0)
        self.ltShell[self.shellnr].addWidget(lblSigma, 2, 0)
        self.ltShell[self.shellnr].addWidget(lblC3, 3, 0)
        self.ltShell[self.shellnr].addWidget(lblC4, 4, 0)
        self.ltShell[self.shellnr].addWidget(lblC5, 5, 0)
        self.ltShell[self.shellnr].addWidget(lblC6, 6, 0)
        self.ltShell[self.shellnr].addWidget(lblE0, 7, 0)  
        self.ltShell[self.shellnr].addWidget(lblAmp, 8, 0) 
        self.ltShell[self.shellnr].addWidget(lblPha, 9, 0) 
        
        for i in range(5):           
            self.ltShell[self.shellnr].addWidget(self.shellN[self.shellnr][i], 0, 2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellR[self.shellnr][i], 1,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellSigma[self.shellnr][i], 2,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC3[self.shellnr][i], 3,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC4[self.shellnr][i], 4,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC5[self.shellnr][i], 5,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellC6[self.shellnr][i], 6,  2*i+1)
            self.ltShell[self.shellnr].addWidget(self.shellE0[self.shellnr][i], 7,  2*i+1)
        self.ltShell[self.shellnr].addWidget(self.shellAmp[self.shellnr], 8, 1, 1, 8)
        self.ltShell[self.shellnr].addWidget(self.shellPha[self.shellnr], 9, 1, 1, 8)
        

        
        for j in range(7):
            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Min. limit"), j, 2)
            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Max. limit"), j, 4)
            self.ltShell[self.shellnr].addWidget(QtGui.QLabel("Accuracy"), j, 6)
        
        self.tabs[self.shellnr].setLayout(self.ltShell[self.shellnr])
        
        self.shellnr = self.shellnr +1
        
    def resetFit(self):
        self.chkFirstFit.setChecked(True)
        
    def savefiledialog_qtgui(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setAcceptMode(1) # save dialog
        dlg.setNameFilters(["All files (*.*)"])
#        dlg.setDirectory(self.currentdir)
        if dlg.exec_():
            flist =  dlg.selectedFiles()
            return flist[0]
        else:
            return ""
        

        