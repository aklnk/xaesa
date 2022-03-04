# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 11:20:43 2018

@author: akali
"""
from init import QTVer

if QTVer == 4:
    from PyQt4 import QtGui, QtCore
    
if QTVer == 5:
    from PyQt5 import QtWidgets as QtGui
    from PyQt5 import QtCore

class xaesa_settings(QtGui.QWidget):
    def __init__(self, parent=None):
        super(xaesa_settings, self).__init__()
        
        # self.rebinSmoothing = 1
        self.rebinE1 = 5
        self.rebinE1E0 = 0.1
        self.rebinE0E3 = 0.02
        
        
        self.initUI()

        
    def initUI(self):
        #################  Rebin group ###################################
        grpRebin = QtGui.QGroupBox("Rebin parameters")
        
        # self.lblRebinSmoothing = QtGui.QLabel("Rebin smoothing factor")
        # self.edtRebinSmoothing = QtGui.QLineEdit(str(self.rebinSmoothing))
        self.lblrebinE1 = QtGui.QLabel("Rebin dE (E < E1)")
        self.edtrebinE1 = QtGui.QLineEdit(str(self.rebinE1))
        self.lblrebinE1E0 = QtGui.QLabel("Rebin dE (E1 < E < E0+50)")
        self.edtrebinE1E0 = QtGui.QLineEdit(str(self.rebinE1E0))
        self.lblrebinE0E3 = QtGui.QLabel("Rebin dk (E0+50 < E < E3)")
        self.edtrebinE0E3 = QtGui.QLineEdit(str(self.rebinE0E3))
        loutRebin = QtGui.QGridLayout()
        # loutRebin.addWidget(self.lblRebinSmoothing, 0, 0)
        # loutRebin.addWidget(self.edtRebinSmoothing, 0, 1) 
        loutRebin.addWidget(self.lblrebinE1,        0, 0)
        loutRebin.addWidget(self.edtrebinE1,        0, 1) 
        loutRebin.addWidget(self.lblrebinE1E0,      1, 0)
        loutRebin.addWidget(self.edtrebinE1E0,      1, 1) 
        loutRebin.addWidget(self.lblrebinE0E3,      2, 0)
        loutRebin.addWidget(self.edtrebinE0E3,      2, 1) 
        grpRebin.setLayout(loutRebin)
        
        ################ Apply button ###############
        self.btnApply = QtGui.QPushButton("Apply")
        self.btnApply.clicked.connect(self.apply)
        
        
        #################Main Layout ###############
        loutMain = QtGui.QGridLayout()
        loutMain.addWidget(grpRebin,        0, 0)
        loutMain.addWidget(self.btnApply,   1, 0)
        
        self.setLayout(loutMain)
        
    def apply(self):
        # self.rebinSmoothing = float(self.edtRebinSmoothing.text())
        self.rebinE1 = float(self.edtrebinE1.text())
        self.rebinE1E0 = float(self.edtrebinE1E0.text())
        self.rebinE0E3 = float(self.edtrebinE0E3.text())
   

#class TestWindow(QtGui.QMainWindow):
#    def __init__(self, parent=None):
#        super(TestWindow, self).__init__()  
#        
#        self.wid = xaesa_settings(self)
#        self.setCentralWidget(self.wid) 
#
#        self.show()
#        
#    def closeEvent(self, event):
#        super(TestWindow, self).closeEvent(event)
#        print('closed')     
#        
#if __name__ == '__main__':
#    app = QtGui.QApplication(argv)
#
#    main = TestWindow()
#
#    exit(app.exec_())