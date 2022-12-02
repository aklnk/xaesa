# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 11:20:43 2018

@author: akali
"""
from init import QTVer

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

class xaesa_viewer(QtGui.QWidget):
    def __init__(self, parent=None):
        super(xaesa_viewer, self).__init__()
        
        self.x_data = [] #list with arrays
        self.y_data = []
        
        self.labels = [] #list with labels for each line
        
        self.x_caption = ""
        self.y_caption = "" 

        # mode = 0 - each line new color  
        # mode = 1 - combine 2 lines with the same color (for FT comparison)
        self.mode = 0      

        self.plotXRange = [] #two element [min, max] (for XANES)    
        self.horizontalLine = []
        
        self.lines = []
        
        self.initUI()

        
    def initUI(self):
        #Figures 
        self.fig = plt.figure(2, figsize=(15, 6))
        # self.ax_exafs = self.fig.add_subplot(111)

        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        
        fnt = self.tbar.font()
        fnt.setPointSize(20)
        self.tbar.setFont(fnt)
        
#        plt.tight_layout()    
        
        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
              
        self.setLayout(lfig)
        
        self.canv.draw()
        
    def plot(self):
        
        plt.clf()
        self.ax_exafs = self.fig.add_subplot(111)
        
        
        if self.mode == 0: #compare exafs        
            for i in range(len(self.x_data)):
                l, = self.ax_exafs.plot(self.x_data[i], self.y_data[i], label = self.labels[i])
                self.ax_exafs.set_xlabel(self.x_caption)
                self.ax_exafs.set_ylabel(self.y_caption)
#                self.ax_exafs.set_xlabel('Wavevector k, $\AA^{-1}$')
#                self.ax_exafs.set_ylabel('EXAFS, $\AA^{-2}$')
                if self.horizontalLine != []:
                    self.ax_exafs.axhline(y=1, linewidth=0.5, color = 'k', linestyle='--',)
                if self.plotXRange != []:
                    self.ax_exafs.set_xlim(self.plotXRange)
                self.lines.append(l)
                
        if self.mode == 1: #compare ft       
            for i in range(len(self.r)/2):
                line1,  = self.ax_exafs.plot(self.x_data[2*i], self.y_data[2*i], label = self.labels[i])
                line2,  = self.ax_exafs.plot(self.x_data[2*i+1], self.y_data[2*i+1])
                line2.set_color(line1.get_color())
                line2.set_linestyle('dotted')
                self.ax_exafs.set_xlabel(self.x_caption)
                self.ax_exafs.set_ylabel(self.y_caption)
#                self.ax_exafs.set_xlabel('Distance R, $\AA$')
#                self.ax_exafs.set_ylabel('Fourier transform, $\AA^{-3}$')
                self.lines.append(line1)
                self.lines1.append(line2)
                      
                
                
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
   

class xaesaViewerWindow(QtGui.QDialog):
    def __init__(self, parent=None):
        super(xaesaViewerWindow, self).__init__()  
        
        self.viewer = xaesa_viewer(self)
        
        self.btnCancel = QtGui.QPushButton('Exit')
        self.btnCancel.clicked.connect(self.cancel)
        
        lout = QtGui.QGridLayout()
        lout.addWidget(self.viewer)
        lout.addWidget(self.btnCancel)
        self.setLayout(lout)
#        self.setCentralWidget(self.wid) 

#        self.show()
        
    def closeEvent(self, event):
        super(xaesaViewerWindow, self).closeEvent(event)    
        
    def cancel(self):
        plt.close(2) 
        self.close()
#        
#if __name__ == '__main__':
#    app = QtGui.QApplication(argv)
#
#    main = TestWindow()
#
#    exit(app.exec_())