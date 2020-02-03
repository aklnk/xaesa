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

from numpy import where, logical_and, sum, transpose, amax, concatenate, savetxt, asarray, array

class xaesa_rxes(QtGui.QWidget):
    def __init__(self, parent=None):
        super(xaesa_rxes, self).__init__()
        
        self.fluoEnergy = []
        self.incidentEnergy = []
        self.rxes = []
        self.energyTransfer = False
        
        self.initUI()

        
    def initUI(self):
        #Figures 
        self.fig = plt.figure(997, figsize=(15, 6))
        self.ax_rxes = plt.subplot2grid((1,2), (0,0))
        self.ax_cuts = plt.subplot2grid((1,2), (0,1))
        
        if self.energyTransfer:
            self.ax_rxes.set_xlabel('Incident energy, eV')
            self.ax_rxes.set_ylabel('Energy transfer, eV')
            self.ax_cuts.set_xlabel('Incident energy, eV')
            self.ax_cuts.set_ylabel('Intensity, a.u.')
            
        else:
            self.ax_rxes.set_xlabel('Incident energy, eV')
            self.ax_rxes.set_ylabel('Emission energy, eV')
            self.ax_cuts.set_xlabel('Incident energy, eV')
            self.ax_cuts.set_ylabel('HERFD XANES, a.u.')
            
        
        self.fig.canvas.mpl_connect('button_press_event', self.imageDoubleClick)

        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        
        fnt = self.tbar.font()
        fnt.setPointSize(20)
        self.tbar.setFont(fnt)
        
#        plt.tight_layout()    
        
        self.sl = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.sl.setMinimum(0)
        self.sl.setMaximum(100)
        self.sl.setValue(100)
        self.sl.setTickPosition(QtGui.QSlider.TicksBelow)
        self.sl.valueChanged.connect(self.valuechange)
#        self.sl.setTickInterval(5)
        
        glout = QtGui.QGridLayout()
        
        glout.addWidget(self.sl, 0, 0, 1, 2)
        
        glout.addWidget(QtGui.QLabel("Energy for HERFD-XANES"), 1, 0 )
        glout.addWidget(QtGui.QLabel("Delta energy for HERFD-XANES"), 2, 0 )
        
        self.edtEnergyHerfd = QtGui.QLineEdit()
        glout.addWidget(self.edtEnergyHerfd , 1, 1 )
        
        self.edtDeltaHerfd = QtGui.QLineEdit()
        glout.addWidget(self.edtDeltaHerfd , 2, 1 )
        
        self.btnCalculateHerfd = QtGui.QPushButton("Show HERFD XANES")
        self.btnCalculateHerfd.clicked.connect(self.showHerfd)
        glout.addWidget(self.btnCalculateHerfd , 3, 0, 1, 2 )
        
        self.btnCalculateHerfd = QtGui.QPushButton("save HERFD XANES")
        self.btnCalculateHerfd.clicked.connect(self.saveHerfd)
        glout.addWidget(self.btnCalculateHerfd , 4, 0, 1, 2 )
        
        glout.addWidget(QtGui.QLabel("Energy transfer min "), 1, 2 )
        glout.addWidget(QtGui.QLabel("Energy transfer max"), 2, 2 )
        
        self.edtEnergyTransMin = QtGui.QLineEdit()
        glout.addWidget(self.edtEnergyTransMin , 1, 3 )
        
        self.edtEnergyTransMax = QtGui.QLineEdit()
        glout.addWidget(self.edtEnergyTransMax , 2, 3 )
        
        self.btnCalculateHerfd = QtGui.QPushButton("Show energy transfer")
        self.btnCalculateHerfd.clicked.connect(self.showEnergyTrans)
        glout.addWidget(self.btnCalculateHerfd , 3, 2, 1, 2 )
        
        self.btnCalculateHerfd = QtGui.QPushButton("Save energy transfer")
        self.btnCalculateHerfd.clicked.connect(self.saveEnergyTrans)
        glout.addWidget(self.btnCalculateHerfd , 4, 2, 1, 2 )
        
        
        
        self.btnCancel = QtGui.QPushButton('Exit')
        self.btnCancel.clicked.connect(self.cancel)
        
        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
        
        lfig.addLayout(glout)
        
        lfig.addWidget(self.btnCancel)
              
        self.setLayout(lfig)
        
        self.canv.draw()
        
    def plot(self):
        
        savetxt("x.dat", self.incidentEnergy)
        savetxt("y.dat", self.fluoEnergy)
        savetxt("z.dat", self.rxes)
        
#        self.image = self.ax_rxes.imshow(self.rxes, vmin = 0, vmax = self.zmax, aspect='auto', interpolation='none',cmap='gist_ncar')
        self.image = self.ax_rxes.pcolormesh(self.incidentEnergy, self.fluoEnergy, self.rxes)
        
#        if  amax(self.rxes) > 1:
#            self.sl.setMaximum(amax(self.rxes))
                
        self.fig.tight_layout()
        
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
        
    def showHerfd(self):
        
        emin = float(self.edtEnergyHerfd.text()) - float(self.edtDeltaHerfd.text())
        emax = float(self.edtEnergyHerfd.text()) + float(self.edtDeltaHerfd.text())
        
        if self.energyTransfer: #plot energy transfer intensity versus incident energy
            indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
            tmp_array = []
            for i in range(len(self.incidentEnergy[0])):
                indexes_where = where( logical_and(self.fluoEnergy[:,i]>emin, self.fluoEnergy[:,i]<emax) )
                print(self.rxes[indexes_where, i][0])
                number = sum(self.rxes[indexes_where, i][0]) / len(self.rxes[indexes_where, i][0])
                print(number)
                tmp_array.append(number)

#            hsum = sum(tmp_array, axis=0)

            self.ax_cuts.clear()
            self.ax_cuts.set_xlabel('Incident energy, eV')
            self.ax_cuts.set_ylabel('Intensity, a.u.')
            self.ax_cuts.plot(self.incidentEnergy[0], tmp_array)

            
        else: # plot normal herfd           
            
            indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
    
            hsum = sum( (self.rxes[indexes_where, :]), axis=1)
    
            self.ax_cuts.clear()
    
            self.ax_cuts.set_xlabel('Incident energy, eV')
            self.ax_cuts.set_ylabel('HERFD XANES, a.u.')
    
            self.ax_cuts.plot(self.incidentEnergy[0], hsum[0])
        
        self.canv.draw()
        
    def saveHerfd(self):
        
        emin = float(self.edtEnergyHerfd.text()) - float(self.edtDeltaHerfd.text())
        emax = float(self.edtEnergyHerfd.text()) + float(self.edtDeltaHerfd.text())
        
        if self.energyTransfer: #plot energy transfer intensity versus incident energy
            indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
            tmp_array = []
            for i in range(len(self.incidentEnergy[0])):
                indexes_where = where( logical_and(self.fluoEnergy[:,i]>emin, self.fluoEnergy[:,i]<emax) )
                print(self.rxes[indexes_where, i][0])
                number = sum(self.rxes[indexes_where, i][0]) / len(self.rxes[indexes_where, i][0])
                print(number)
                tmp_array.append(number)

#            hsum = sum(tmp_array, axis=0)

            forSave = concatenate(([self.incidentEnergy[0]], [tmp_array]))           

            
        else: # plot normal herfd
        
            indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
            
    
            hsum = sum( (self.rxes[indexes_where, :]), axis=1)
            
            forSave = concatenate(([self.incidentEnergy[0]], [hsum[0]]))
            
        fn = self.savefiledialog_qtgui()
        if fn == "":
            return
        
        savetxt(fn, transpose(forSave) )
        
    def showEnergyTrans(self):
        
        emin = float(self.edtEnergyTransMin.text())
        emax = float(self.edtEnergyTransMax.text())
        
        tmp_array = []
        #vert sum Energy transfer
        indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
        nr_points = len(indexes_where[0])
        for i in range(len(self.incidentEnergy[0])):
            indexes_where = where( logical_and(self.fluoEnergy[:,i]>emin, self.fluoEnergy[:,i]<emax) )
            first_i = indexes_where[0][0]
            i_to_take = array(range(first_i, first_i+nr_points))
            tmp_array.append(self.rxes[[i_to_take], i][0])
            
            
        hsum = sum(tmp_array, axis=0)

        indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )

        self.ax_cuts.clear()
        self.ax_cuts.set_xlabel('Incident energy, eV')
        self.ax_cuts.set_ylabel('Intensity, a.u.')
        self.ax_cuts.plot(self.fluoEnergy[:,0][indexes_where], hsum)
        self.canv.draw()
        
    def saveEnergyTrans(self):
        
        emin = float(self.edtEnergyTransMin.text())
        emax = float(self.edtEnergyTransMax.text())
        
        tmp_array = []
        #vert sum Energy transfer
        indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
        nr_points = len(indexes_where[0])
        for i in range(len(self.incidentEnergy[0])):
            indexes_where = where( logical_and(self.fluoEnergy[:,i]>emin, self.fluoEnergy[:,i]<emax) )
            first_i = indexes_where[0][0]
            i_to_take = array(range(first_i, first_i+nr_points))
            tmp_array.append(self.rxes[[i_to_take], i][0])
        hsum = sum(tmp_array, axis=0)

        indexes_where = where( logical_and(self.fluoEnergy[:,0]>emin, self.fluoEnergy[:,0]<emax) )
        
        forSave = concatenate(([self.fluoEnergy[:,0][indexes_where]], [hsum]))
        
        fn = self.savefiledialog_qtgui()
        if fn == "":
            return
        
        savetxt(fn, transpose(forSave) )

        
    def valuechange(self):
#        size = self.sl.value()
        print(self.sl.value())
        
        self.image.set_clim(vmin = 0, vmax = (self.sl.value()/100)*amax(self.rxes) )
        self.canv.draw()
        
    def imageDoubleClick(self, event):
        print(event)
        if event.dblclick == True:
            #ask to change scale
            num, ok = QtGui.QInputDialog.getDouble(self,
                                               "Z scale",
                                               "Z max",
                                               value = self.zmax, 
                                               max = 20000000,
                                               min = 0,
                                               decimals = 4)
            if ok:
                self.zmax = num
                self.image.set_clim(vmax=self.zmax)
                self.canv.draw()
            else:
                return
            
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
        
        
   

class xaesaRxesWindow(QtGui.QDialog):
    def __init__(self, parent=None):
        super(xaesaRxesWindow, self).__init__()  
        
        self.viewer = xaesa_rxes(self)
        lout = QtGui.QGridLayout()
        lout.addWidget(self.viewer)
        self.setLayout(lout)
#        self.setCentralWidget(self.wid) 

#        self.show()
        
    def closeEvent(self, event):
        super(xaesaRxesWindow, self).closeEvent(event)    
#        
#if __name__ == '__main__':
#    app = QtGui.QApplication(argv)
#
#    main = TestWindow()
#
#    exit(app.exec_())