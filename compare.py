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


class CompareWindow(QtGui.QDialog):

    def __init__(self):
        super(CompareWindow, self).__init__()
        
        self.exafs = []
        self.k = []

        self.r = []
        self.fr = []
        self.fi = []

        self.bftk = []
        self.bftexafs = []

        self.labels = []

        self.energy = []
        self.mju = []
        self.xes = []

        self.E0 = 0
        
        self.lines = []
        self.lines1 = []
        #mode
        #0 - exafs
        #1 - ft
        #2 - bft
        #3 - mju
        #10 - Xes original
        #11 - XES area Normalized
        #12 - XES Max Normalized
        self.mode = 0

        self.initUI()

    def initUI(self):
        
        #Figures 
        self.fig = plt.figure(1, figsize=(15, 6))
        # self.ax_exafs = self.fig.add_subplot(111)

        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        
        fnt = self.tbar.font()
        fnt.setPointSize(20)
        self.tbar.setFont(fnt)
        
        # plt.tight_layout()    
        
        self.btnCancel = QtGui.QPushButton('Exit')
        self.btnCancel.clicked.connect(self.cancel)
        
        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
        
        lfig.addWidget(self.btnCancel)
              
        self.setLayout(lfig)
        
        self.canv.draw()
        
        #wid.setLayout(lfig)
        
    def plot(self):
        
        # self.ax_exafs.clear()
        plt.clf()
        self.ax_exafs = self.fig.add_subplot(111)
        
        
        
        if self.mode == 0: #compare exafs        
            for i in range(len(self.k)):
                l, = self.ax_exafs.plot(self.k[i], self.exafs[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Wavevector k, $\AA^{-1}$')
                self.ax_exafs.set_ylabel('EXAFS, $\AA^{-2}$')
                self.lines.append(l)
                
        if self.mode == 1: #compare ft       
            for i in range(len(self.r)):
                line1,  = self.ax_exafs.plot(self.r[i], self.fr[i], label = self.labels[i], linewidth=1)
                line2,  = self.ax_exafs.plot(self.r[i], self.fi[i], linewidth=1)
                line2.set_color(line1.get_color())
                line2.set_linestyle('dotted')
                self.ax_exafs.set_xlabel('Distance R, $\AA$')
                self.ax_exafs.set_ylabel('Fourier transform, $\AA^{-3}$')
                self.lines.append(line1)
                self.lines1.append(line2)
                
                
        if self.mode == 2: #compare bft       
            for i in range(len(self.bftk)):
                l, = self.ax_exafs.plot(self.bftk[i], self.bftexafs[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Wavevector k, $\AA^{-1}$')
                self.ax_exafs.set_ylabel('EXAFS, $\AA^{-2}$')
                self.lines.append(l)
                
        if self.mode == 3: #compare mju      
            for i in range(len(self.energy)):
                l, = self.ax_exafs.plot(self.energy[i], self.mju[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Energy, eV')
                self.ax_exafs.set_ylabel('Absorption, a.u.')
                self.lines.append(l)
                
                
        if self.mode == 4: #compare xanes     
            for i in range(len(self.energy)):
                l, = self.ax_exafs.plot(self.energy[i], self.mju[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Energy, eV')
                self.ax_exafs.set_ylabel('Absorption, a.u.')
                self.ax_exafs.axhline(y=1, linewidth=0.5, color = 'k', linestyle='--',)
                self.ax_exafs.set_xlim([self.E0-75,self.E0+200])
                self.lines.append(l)
                
        if self.mode == 10: # XES original     
            for i in range(len(self.energy)):
                l, = self.ax_exafs.plot(self.energy[i], self.xes[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Energy, eV')
                self.ax_exafs.set_ylabel('Intensity, a.u.')
                self.lines.append(l)
                
        if self.mode == 11: # XES area normalized
            for i in range(len(self.energy)):
                l, = self.ax_exafs.plot(self.energy[i], self.xes[i], label = self.labels[i], linewidth=1)
                self.ax_exafs.set_xlabel('Energy, eV')
                self.ax_exafs.set_ylabel('Area normalized intensity')
                self.lines.append(l)
                
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
                
        self.fig.tight_layout()
                
        box = self.ax_exafs.get_position()
        self.ax_exafs.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        
        leg = self.ax_exafs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        # we will set up a dict mapping legend line to orig line, and enable
        # picking on the legend line
        self.lined = dict()
        for legline, origline in zip(leg.get_lines(), self.lines):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined[legline] = origline

        if self.mode == 1:
            self.lined1 = dict()
            for legline, origline in zip(leg.get_lines(), self.lines1):
#                legline.set_picker(5)  # 5 pts tolerance
                self.lined1[legline] = origline
        
        self.canv.draw()

        
    def cancel(self):
        #do whatever you need with self.roiGroups  
        self.close()
        
    def onpick(self, event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = self.lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        if self.mode == 1:
            origline = self.lined1[legline]
            origline.set_visible(vis)
        
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        self.fig.canvas.draw()
        
                
        
        