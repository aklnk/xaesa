# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:07:03 2016

@author: sasha
"""

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

from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

from . import xaesa_ft
from .xaesa_ft import FT, BFT, BFTWindow, GETPHASE

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


class DGWindow(QtGui.QDialog):

    def __init__(self):
        super(DGWindow, self).__init__()
        
        self.exafs = []
        self.k = []

        self.exafsdg = [1,2,3]

        self.initUI()

    def initUI(self):
        
        #3wid = QtGui.QWidget(self)
        #self.setCentralWidget(wid)
        
        self.lblx1 = QtGui.QLabel("X1")
        self.lblx2 = QtGui.QLabel("X2")
        self.lblx3 = QtGui.QLabel("x3")
        self.lblx4 = QtGui.QLabel("X4")

        self.edtx1 = QtGui.QLineEdit("")
        self.edtx2 = QtGui.QLineEdit("")
        self.edtx3 = QtGui.QLineEdit("")
        self.edtx4 = QtGui.QLineEdit("")


        lregions = QtGui.QGridLayout()
        lregions.addWidget(self.lblx1, 0, 0)
        lregions.addWidget(self.lblx2, 0, 1)
        lregions.addWidget(self.lblx3, 0, 2)
        lregions.addWidget(self.lblx4, 0, 3)


        lregions.addWidget(self.edtx1, 1, 0)
        lregions.addWidget(self.edtx2, 1, 1)
        lregions.addWidget(self.edtx3, 1, 2)
        lregions.addWidget(self.edtx4, 1, 3)

        self.lblGlitchCriteria = QtGui.QLabel("Glitch criteria")
        self.lblrmax = QtGui.QLabel("R max")
        
        self.edtGlitchCriteria = QtGui.QLineEdit("0.05")
        self.edtrmax = QtGui.QLineEdit("10")
        
        self.lstGlitchList = QtGui.QListWidget()
        self.lstGlitchList.itemClicked.connect(self.lstGlitchListItemClicked)
        
        lregions.addWidget(self.lblGlitchCriteria, 2, 0)
        lregions.addWidget(self.edtGlitchCriteria, 2, 1)
        lregions.addWidget(self.lblrmax, 2, 2)
        lregions.addWidget(self.edtrmax, 2, 3)
        
        lregions.addWidget(self.lstGlitchList, 3, 0, 1, 4)

        
        #Figures 
        self.fig = plt.figure(2, figsize=(12, 6))
        self.ax_exafs = plt.subplot2grid((1,2), (0,1))
        self.ax_exafsdg = plt.subplot2grid((1,2), (0,0))
        self.ax_exafsdg.set_title('Left mouse button to select region before glitch\nRight mouse button to select region after glitch')

        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        
        # set useblit True on gtkagg for enhanced performance
        self.span = SpanSelector(self.ax_exafsdg, self.onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'), button = 1, span_stays=True)
        
        self.span1= SpanSelector(self.ax_exafsdg, self.onselect1, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='green'), button = 3, span_stays=True)

        self.exafs_line, = self.ax_exafs.plot([], [])
#        self.bftexafs_line, = self.ax_exafs.plot([], [])
        self.difference_line, = self.ax_exafs.plot([], [], "o")
        self.exafsdgcomp_line, = self.ax_exafs.plot([], [])
        
        self.exafsdg_line, = self.ax_exafsdg.plot([], [])
        self.difference_line.set_markersize(2)
        
        self.exafsdgcomp_line.set_color('r')
        
        self.glitch_lines = []
        self.minmaxlines = []
        
        self.fig.tight_layout()
        
        self.btnGlitchFinder = QtGui.QPushButton('Glitch finder')
        self.btnGlitchFinder.clicked.connect(self.GlitchFinder)
        
        self.btnDG = QtGui.QPushButton('Deglitch')
        self.btnDG.clicked.connect(self.DG)
        
        self.btnApply = QtGui.QPushButton('Apply')
        self.btnApply.clicked.connect(self.apply)
        
        self.btnCancel = QtGui.QPushButton('Cancel')
        self.btnCancel.clicked.connect(self.cancel)
        
        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
        
        lfig.addLayout(lregions)
        lfig.addWidget(self.btnGlitchFinder)
        lfig.addWidget(self.btnDG)
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
        self.exafs_line.set_xdata(self.k)
        self.exafs_line.set_ydata(self.exafs)
        self.exafsdg_line.set_xdata(self.k)
        self.exafsdg_line.set_ydata(self.exafsdg)
        self.exafsdgcomp_line.set_xdata(self.k)
        self.exafsdgcomp_line.set_ydata(self.exafsdg)
        self.ax_exafs.relim()
        self.ax_exafs.autoscale()
        self.ax_exafsdg.relim()
        self.ax_exafsdg.autoscale()
        self.canv.draw()
        
    def apply(self):
        #do whatever you need with self.roiGroups    
        self.accept()
        
    def cancel(self):
        #do whatever you need with self.roiGroups    
        self.close()
        
    def DG(self):
        if self.edtx1.text() == "" or self.edtx2.text() == ""  \
                          or self.edtx3.text() == "" or self.edtx4.text() == "" :
                              
            return
            
        x1 = float(self.edtx1.text())
        x2 = float(self.edtx2.text())
        x3 = float(self.edtx3.text())
        x4 = float(self.edtx4.text())
        
        
        k_dg = np.array([])
        exafs_dg = np.array([])
        
        for i in range( len(self.k)):
            if( (self.k[i] > x1 and self.k[i] < x2) or (self.k[i] > x3 and self.k[i] < x4) ) :
                k_dg = np.append(k_dg, self.k[i])
                exafs_dg = np.append(exafs_dg, self.exafs[i])
                
        spl = UnivariateSpline(k_dg, exafs_dg)
        
        for i in range( len(self.k)):
            if (self.k[i] > x2 and self.k[i] < x3): #change glitch to spline values
                 self.exafsdg[i] = spl(self.k[i])
        self.plot()
        
    def GlitchFinder(self):
        window = Window_Gauss10(self.k, self.k[0], self.k[len(self.k)-1])

        wexafs = self.exafsdg * window
        
        rmax = float(self.edtrmax.text())

        r, fr, fi = FT(self.k, wexafs, 0, rmax, 0.02)

        bftw = BFTWindow(r, 0, rmax, 0.1)
        forbftim = fi * bftw
        forbftre = fr * bftw
        bftk,  bftefr, bftefi =  BFT(r, forbftre, forbftim, self.k[0], self.k[len(self.k)-1], 0.05)

        wind =  Window_Gauss10(bftk, self.k[0], self.k[len(self.k)-1])
        bftamp = np.sqrt(bftefr*bftefr / (wind*wind) + bftefi*bftefi / (wind*wind))

        bftpha = GETPHASE(bftefr, bftefi)

        bftexafs = bftamp * np.sin(bftpha)
        
        spl = InterpolatedUnivariateSpline(bftk, bftexafs)
        bftexafs_k = spl(self.k)
        
        difference = self.exafsdg - bftexafs_k
        difference2 = difference * difference
        difference = np.sqrt(difference2)
        
        criteria = float(self.edtGlitchCriteria.text())
        
        where_criteria = np.where(difference > criteria)
#        k_criteria = self.k[where_criteria]
        #k_criteria = np.delete(k_criteria, len(k_criteria)-1)

#        self.ax_exafs.clear()
#        self.span = SpanSelector(self.ax_exafs, self.onselect, 'horizontal', useblit=True,
#                    rectprops=dict(alpha=0.5, facecolor='red'), button = 1, span_stays=True)
#        
#        self.span1= SpanSelector(self.ax_exafs, self.onselect1, 'horizontal', useblit=True,
#                    rectprops=dict(alpha=0.5, facecolor='green'), button = 3, span_stays=True)
#        plot1, = self.ax_exafs.plot(self.k, self.exafs)
        glitch_pos = []
#        for i in range(len(where_criteria[0])-2):
#            if where_criteria[0][i+1] == where_criteria[0][i]+1:
#                glitch_pos.append( (k_criteria[i]+k_criteria[i+1])/2)
#                self.ax_exafs.axvline((k_criteria[i]+k_criteria[i+1])/2, color='r', linestyle='--', lw=1)
                
        self.glitch_points = np.split(where_criteria[0], np.where(np.diff(where_criteria[0]) != 1)[0]+1)
        
#        if glitch_points[len(glitch_points)-1] >= self.k(len(self.k)-1):
#            glitch_points.delete(len(glitch_points)-1)
        
        for i in range(len(self.glitch_lines)):
            self.glitch_lines[i].remove()
        
        self.glitch_lines = []
        remove1st = False
        removelast = False
        for i in range(len(self.glitch_points)):
            val = np.sum( self.k[self.glitch_points[i]]) / len(self.glitch_points[i])
            if val+0.5 >= self.k[len(self.k)-1]:
                removelast = True
                continue
            if val-0.5 <= self.k[0]:
                remove1st = True
                continue
            glitch_pos.append( val )
            line = self.ax_exafsdg.axvline(val, color='r', linestyle='--', lw=1)
            self.glitch_lines.append(line)
        
        if remove1st:
            del self.glitch_points[0]
        if removelast:
            del self.glitch_points[len(self.glitch_points)-1]
        
        
        self.lstGlitchList.clear()
        self.lstGlitchList.addItems(list(map(str,glitch_pos))) 
        
#        deriv_exafs = np.gradient(self.exafsdg)
        
        
#        self.bftexafs_line.set_xdata(self.k)
#        self.bftexafs_line.set_ydata(bftexafs_k)
        self.difference_line.set_xdata(self.k)
        self.difference_line.set_ydata(difference)
        
        self.ax_exafs.relim()
        self.ax_exafs.autoscale()
        self.ax_exafsdg.relim()
        self.ax_exafsdg.autoscale()
        self.canv.draw()

        
    def lstGlitchListItemClicked(self):
        
        current =  self.lstGlitchList.currentRow()
        xmax1 = self.k[self.glitch_points[current][0]]-0.05
        xmin1 = xmax1 - 0.2
        
        self.edtx1.setText(str(xmin1))
        self.edtx2.setText(str(xmax1))
        
        
        xmin2 = self.k[self.glitch_points[current][len(self.glitch_points[current])-1]]+0.05
        xmax2 = xmin2 + 0.2
        self.edtx3.setText(str(xmin2))
        self.edtx4.setText(str(xmax2))
        
        for i in range(len(self.minmaxlines)):
            self.minmaxlines[i].remove()
            
        self.minmaxlines = []
        
        line = self.ax_exafsdg.axvline(xmin1, color='g', linestyle='--', lw=1)
        self.minmaxlines.append(line)
        line = self.ax_exafsdg.axvline(xmax1, color='g', linestyle='--', lw=1)
        self.minmaxlines.append(line)
        line = self.ax_exafsdg.axvline(xmin2, color='g', linestyle='--', lw=1)
        self.minmaxlines.append(line)
        line = self.ax_exafsdg.axvline(xmax2, color='g', linestyle='--', lw=1)
        self.minmaxlines.append(line)
        self.canv.draw()
        
        
        
        
        
        