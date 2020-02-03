# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:07:03 2016

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

from scipy import linalg as LA

from matplotlib.widgets import SpanSelector

from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

#from . import ft
#from . import xaslib

#from tkinter import messagebox

from scipy.optimize import least_squares

from sklearn.decomposition import PCA

#from tkinter.filedialog import  asksaveasfilename
#from tkinter import Tk

def lin_comb_func(x, basis, fit):
    
#    print(basis * x[:,None])
#    print(np.sum(basis * x[:,None], axis=0))
    
    return np.sum(basis * x[:,None], axis=0) - fit

class PCAWindow(QtGui.QDialog):

    def __init__(self):
        super(PCAWindow, self).__init__()
        
        self.mainform = None
        
        self.mju_exafs = 1 # 0 - mju, 1 - exafs, 2 - xes

#        self.initUI()

    def initUI(self):
        
        self.fit_result = []
        self.fit_basis = []

        lblBasis = QtGui.QLabel("Select basis:")
        lblFit = QtGui.QLabel("Select spectra to fit (multiple selection possible):")
        lblResult = QtGui.QLabel("Fit coeficients:")
        lblNComponents = QtGui.QLabel("Max number of PCA components:")
        
        self.edtNComponents = QtGui.QLineEdit("3")
        
        self.lstBasis = QtGui.QListWidget()
        self.lstFit = QtGui.QListWidget()
        self.lstBasis.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.lstFit.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.lstFit.itemClicked.connect(self.lstFitItemClicked)
        
        self.lstFitResult = QtGui.QListWidget()
        
        anfn = []
        fit_res = []
        for i in range(0, self.mainform.lstSpectra.count()):
            anfn.append(str(self.mainform.lstSpectra.item(i).text()))
            fit_res.append("")
            self.fit_result.append([])
            self.fit_basis.append([])
            
        self.lstBasis.addItems(anfn)
        self.lstFit.addItems(anfn)
        self.lstFitResult.addItems(fit_res)

        lout = QtGui.QGridLayout()
        lout.addWidget(lblBasis, 0,0)
        lout.addWidget(lblFit,0,1 )
        lout.addWidget(lblResult,0,2)
        lout.addWidget(self.lstBasis, 1, 0)
        lout.addWidget(self.lstFit, 1, 1)
        lout.addWidget(self.lstFitResult,1,2)
        lout.addWidget(lblNComponents, 2, 0)
        lout.addWidget(self.edtNComponents, 2, 1)
        
       
        #Figures 
        self.fig = plt.figure(6, figsize=(12, 6))
        self.ax_exafs = plt.subplot2grid((1,2), (0,1))
        self.ax_fit = plt.subplot2grid((1,2), (0,0))

        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        

        self.exafs_line, = self.ax_exafs.plot([], [])
#        self.bftexafs_line, = self.ax_exafs.plot([], [])
        self.difference_line, = self.ax_exafs.plot([], [], "o")
        self.exafsdgcomp_line, = self.ax_exafs.plot([], [])
        
#        self.fit_line, = self.ax_fit.plot([], [], "o")
#        self.fit_line.set_markersize(2)
        
        self.exafsdgcomp_line.set_color('r')
        
        self.glitch_lines = []
        self.minmaxlines = []
        
        self.fig.tight_layout()
        
        self.btnFit = QtGui.QPushButton('Find PCA components')
        self.btnFit.clicked.connect(self.findPcaComponents)
        lout.addWidget(self.btnFit, 2, 2)
        
        self.btnPCATransform = QtGui.QPushButton('PCA Transform')
        self.btnPCATransform.clicked.connect(self.PCATransform)
        
        self.btnApply = QtGui.QPushButton('Save fit results...')
        self.btnApply.clicked.connect(self.save_fit)
        
        self.btnCancel = QtGui.QPushButton('Exit')
        self.btnCancel.clicked.connect(self.cancel)
        
        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
        
        lfig.addLayout(lout)
        lfig.addWidget(self.btnPCATransform)
        lfig.addWidget(self.btnApply)
        lfig.addWidget(self.btnCancel)
              
        self.setLayout(lfig)
        
        #wid.setLayout(lfig)
        
    def onselect(self, xmin, xmax):
        self.edtx1.setText(str(xmin))
        self.edtx2.setText(str(xmax))
        
    def onselect1(self, xmin, xmax):
        self.edtx3.setText(str(xmin))
        self.edtx4.setText(str(xmax))

        
    def plot(self):
#        self.ax_exafs.clear()
#        # set useblit True on gtkagg for enhanced performance
#        self.span = SpanSelector(self.ax_exafs, self.onselect, 'horizontal', useblit=True,
#                    rectprops=dict(alpha=0.5, facecolor='red'), button = 1, span_stays=True)
#        
#        self.span1= SpanSelector(self.ax_exafs, self.onselect1, 'horizontal', useblit=True,
#                    rectprops=dict(alpha=0.5, facecolor='green'), button = 3, span_stays=True)
#        self.ax_exafsdg.clear()
#        self.ax_exafs.plot(self.k, self.exafs)
#        self.ax_exafsdg.plot(self.k, self.exafsdg)
#        self.exafs_line.set_xdata(self.k)
#        self.exafs_line.set_ydata(self.exafs)
#        self.exafsdg_line.set_xdata(self.k)
#        self.exafsdg_line.set_ydata(self.exafsdg)
#        self.exafsdgcomp_line.set_xdata(self.k)
#        self.exafsdgcomp_line.set_ydata(self.exafsdg)
#        self.ax_exafs.relim()
#        self.ax_exafs.autoscale()
#        self.ax_exafsdg.relim()
#        self.ax_exafsdg.autoscale()
        self.canv.draw()
        
    def save_fit(self):

        header = "#"
        
        
        sa = [[],[]]
        fmt="%10s"
        for i in range(self.lstFit.count()):
            if self.fit_basis[i] != []:
                sa[0].append([self.lstFit.item(i).text()])
                sa[1].append(self.fit_result[i])
#                for i in range(len(self.fit_basis[i])):
#                    fmt = fmt + " %10.3f"

        
#        Tk().withdraw()
#        fn = asksaveasfilename()
#        if fn == "":
#            return

        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setAcceptMode(1) # save dialog
        dlg.setNameFilters(["All files (*.*)"])
        dlg.setDirectory(os.getcwd())
        if dlg.exec_():
            flist =  dlg.selectedFiles()
            fn =  flist[0]
        else:
            return

        print(sa)
        print(fmt)
        print(np.hstack((sa[0], sa[1])))
        sa[0] = np.array(sa[0])
        sa[1] = np.array(sa[1])
#        sa[1] = reshape(1, len(sa[0])*len(self.fit_basis[0]))
        np.savetxt(fn, np.hstack((sa[0], sa[1])), fmt='%s')
        
    def cancel(self):
   
        self.close()
        

    def findPcaComponents(self):
        sel_basis = self.lstBasis.selectedItems()
        sel_fit = self.lstFit.selectedItems()
        sel_basisi = self.lstBasis.selectedIndexes()
        sel_fiti = self.lstFit.selectedIndexes()
        
        if len(sel_basis) < 1:
#            messagebox.showinfo("Select basis spectra", "Please select at least 1 spectra for basis")
            return
            
#        if len(sel_fit) < 1:
##            messagebox.showinfo("Select fit spectra", "Please select at least 1 spectrum to fit")
#            return
            
        basis = []
        basisi = []
        variables = []
        varbounds = [[],[]]
        for i in range(len(sel_basis)):
            y= sel_basisi[i].row()
            if self.mju_exafs == 3: 
                basis.append(self.mainform.dataClasses[y].xesAreaNorm)
            if self.mju_exafs == 2:
                idx = np.argmin(abs(self.mainform.dataClasses[y].energy - \
                                    self.mainform.dataClasses[y].normalizationEnergy))                            
                basis.append(self.mainform.dataClasses[y].mjuMinusVictoreen / self.mainform.dataClasses[y].mju0[idx])
            if self.mju_exafs == 1: 
                basis.append(self.mainform.dataClasses[y].exafsZLC)
            if self.mju_exafs == 0:
                basis.append(self.mainform.dataClasses[y].mjuMinusVictoreen)
            basisi.append(y)
            variables.append(1)
            varbounds[0].append(0)
            varbounds[1].append(1E6)
            
#        basis = np.transpose(basis)
            
        ncomp = int(self.edtNComponents.text())

        self.pca = PCA(n_components=ncomp)
        self.pca.fit(basis)
        pca_score = self.pca.explained_variance_ratio_
        V = self.pca.components_

#        print(V)
        print("Component score", pca_score)
        
#        tr = self.pca.transform([V[1]*pca_score[1]])
#        print("transform", tr)
        
#        ft = pca.transform([V[1] * pca_score[1]])
#        print(ft)
#        
#        a = 3 * V.T
#        print(a)
        
        #PCA with numpy
#        x = np.transpose(basis)
#
#        #centering the data
#        x -= np.mean(x, axis = 0)  
#        
#        cov = np.cov(x, rowvar = False)
#        
#        evals , evecs = LA.eigh(cov)
#        
#        idx = np.argsort(evals)[::-1]
#        evecs = evecs[:,idx]
#        evals = evals[idx]
#        
#        print(evals)
#        print(evecs)
        
        
#        self.ax_fit.clear()
#        for i in range(len(sel_fit)):
#            y= sel_fiti[i].row()
#            if self.mju_exafs == 2: 
#                fit = self.mainform.dataClasses[y].xesAreaNorm
#            if self.mju_exafs == 1:
#                fit = self.mainform.dataClasses[y].exafsZLC
#            if self.mju_exafs == 0:
#                fit = self.mainform.dataClasses[y].mjuMinusVictoreen
#        
##        lin_comb_func(np.array([1,1]), np.array(basis), np.array(fit))    
#            lsq_result = least_squares(lin_comb_func, np.array(variables), \
##                                   method = 'dogbox',
#                           ftol = 1e-12,  
#                           max_nfev = 1000,
#                           bounds = varbounds,
##                                   tr_solver = 'lsmr',
##                                   jac = '3-point',
##                                   loss='soft_l1',
##                                   f_scale=0.1,
#                            verbose = 2,
##                                    x_scale = [1,1,0.001],
#                           args=(basis, fit) )
#            
#            print (lsq_result.x)
#            self.lstFitResult.item(y).setText(str(lsq_result.x))
#            
#            func = lin_comb_func(lsq_result.x, np.array(basis), np.array(fit))
#
#            self.fit_basis[y] = basisi
#            self.fit_result[y] = lsq_result.x
#
#            if self.mju_exafs == 2:
#                k_energy = self.mainform.dataClasses[y].energy
#            if self.mju_exafs == 1:
#                k_energy = self.mainform.dataClasses[y].k
#            if self.mju_exafs == 0:
#                k_energy = self.mainform.dataClasses[y].energy
#
#            
#            self.exafs_line.set_xdata(k_energy)
#            self.exafs_line.set_ydata(fit)
#            self.exafsdgcomp_line.set_xdata(k_energy)
#            self.exafsdgcomp_line.set_ydata(func+fit)
#            self.ax_exafs.relim()
#            self.ax_exafs.autoscale()
#            
#            self.ax_fit.plot(np.zeros(len(lsq_result.x))+y, lsq_result.x, "o")

        if self.mju_exafs == 3:
            k_energy = self.mainform.dataClasses[y].energy
        if self.mju_exafs == 2:
                k_energy = self.mainform.dataClasses[y].energy
        if self.mju_exafs == 1:
            k_energy = self.mainform.dataClasses[y].k
        if self.mju_exafs == 0:
            k_energy = self.mainform.dataClasses[y].energy  
            
        self.ax_exafs.clear()    
#        for i in range(2):
#            self.ax_exafs.plot(k_energy, V[i]*pca_score[i])
#        sp = self.pca.inverse_transform(tr)
#        print(sp)
#        for i in range(len(sp)):
#            self.ax_exafs.plot(k_energy, sp[i])
        for i in range(len(V)):
            self.ax_exafs.plot(k_energy, V[i]*pca_score[i])
#            
#        self.ax_exafs.plot(k_energy, sp[0])   
#        self.ax_exafs.plot(k_energy, self.mainform.dataClasses[5].exafsZLC)  
#        self.ax_exafs.plot(k_energy, evecs[0])        
        self.ax_exafs.relim()
        self.ax_exafs.autoscale()     
        self.canv.draw()   

    def PCATransform(self):
        self.ax_fit.clear()
        sel_fit = self.lstFit.selectedItems()
        sel_fiti = self.lstFit.selectedIndexes()
        for i in range(len(sel_fit)):
            y= sel_fiti[i].row()
            if self.mju_exafs == 3: 
                fit = self.mainform.dataClasses[y].xesAreaNorm
            if self.mju_exafs == 2:
                idx = np.argmin(abs(self.mainform.dataClasses[y].energy - \
                                    self.mainform.dataClasses[y].normalizationEnergy))                            
                fit = self.mainform.dataClasses[y].mjuMinusVictoreen / self.mainform.dataClasses[y].mju0[idx]
            if self.mju_exafs == 1:
                fit = self.mainform.dataClasses[y].exafsZLC
            if self.mju_exafs == 0:
                fit = self.mainform.dataClasses[y].mjuMinusVictoreen
                
            pca_tr = self.pca.transform([fit]) 
            pca_inv = self.pca.inverse_transform(pca_tr)
                
            self.lstFitResult.item(y).setText(str(pca_tr[0]))
            
            if self.mju_exafs == 3:
                k_energy = self.mainform.dataClasses[y].energy
            if self.mju_exafs == 2:
                k_energy = self.mainform.dataClasses[y].energy
            if self.mju_exafs == 1:
                k_energy = self.mainform.dataClasses[y].k
            if self.mju_exafs == 0:
                k_energy = self.mainform.dataClasses[y].energy
                
#            print(self.pca.score_samples([fit]))

            self.ax_exafs.clear() 
            self.ax_exafs.plot(k_energy, fit)
            self.ax_exafs.plot(k_energy, pca_inv[0])
            self.ax_exafs.relim()
            self.ax_exafs.autoscale()
            self.canv.draw()
         
            
            
    def lstFitItemClicked(self):
        current = self.lstFit.currentRow()
        if self.fit_basis[current] != []:
            for i in range(self.lstBasis.count()):
#                self.lstBasis.setSelected(self.lstBasis.item(i), False)
                self.lstBasis.item(i).setSelected(False)
            for i in self.fit_basis[current]:            
#                self.lstBasis.setSelected(self.lstBasis.item(self.fit_basis[current][i]), True)
                self.lstBasis.item(i).setSelected(True)
                
            basis = []
            basisi = []
            for i in range(len(self.fit_basis[current])):
                if self.mju_exafs == 2:
                    basis.append(self.mainform.dataClasses[self.fit_basis[current][i]].xesAreaNorm)
                if self.mju_exafs == 1:
                    basis.append(self.mainform.dataClasses[self.fit_basis[current][i]].exafsZLC)
                if self.mju_exafs == 0:
                    basis.append(self.mainform.dataClasses[self.fit_basis[current][i]].mjuMinusVictoreen)
                basisi.append(current)
             
            if self.mju_exafs == 2:
                k_energy = self.mainform.dataClasses[current].energy
                abs_exafs = self.mainform.dataClasses[current].xesAreaNorm
            if self.mju_exafs == 1:
                k_energy = self.mainform.dataClasses[current].k
                abs_exafs = self.mainform.dataClasses[current].exafsZLC
            if self.mju_exafs == 0:
                k_energy = self.mainform.dataClasses[current].energy
                abs_exafs = self.mainform.dataClasses[current].mjuMinusVictoreen                
            
            func = lin_comb_func(self.fit_result[current], np.array(basis), np.array(abs_exafs))
            
            
            self.exafs_line.set_xdata(k_energy)
            self.exafs_line.set_ydata(abs_exafs)
            self.exafsdgcomp_line.set_xdata(k_energy)
            self.exafsdgcomp_line.set_ydata(func+abs_exafs)
            self.ax_exafs.relim()
            self.ax_exafs.autoscale()
            self.canv.draw()
            
        else:
            self.exafs_line.set_xdata([])
            self.exafs_line.set_ydata([])
            self.exafsdgcomp_line.set_xdata([])
            self.exafsdgcomp_line.set_ydata([])
            self.ax_exafs.relim()
            self.ax_exafs.autoscale()
            self.canv.draw()
        
      

        
        
        
        