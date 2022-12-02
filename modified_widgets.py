#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 15:02:40 2022

@author: kalinko
"""
from init import QTVer
if QTVer == 4:
    from PyQt4 import QtGui, QtCore

    
if QTVer == 5:
    from PyQt5 import QtWidgets as QtGui
    from PyQt5 import QtCore
    
class QLineEditScroll(QtGui.QLineEdit):
    KEY_shift = QtCore.Qt.Key_Shift
    KEY_ctrl = QtCore.Qt.Key_Control
    
    def __init__(self, base_change_constant = 1, *args, **kwargs):
        QtGui.QLineEdit.__init__(self, *args, **kwargs)
        self.shiftKeyPressed = False
        self.ctrlKeyPressed = False
        self.base_change_constant = base_change_constant
        
    def keyPressEvent(self, event):
        print(event.key())
        if event.key() == QLineEditScroll.KEY_shift:
            self.shiftKeyPressed = True
        if event.key() == QLineEditScroll.KEY_ctrl:
            self.ctrlKeyPressed = True
        QtGui.QLineEdit.keyPressEvent(self, event)

    def keyReleaseEvent(self, event):
        if event.key() == QLineEditScroll.KEY_shift:
            self.shiftKeyPressed = False
        if event.key() == QLineEditScroll.KEY_ctrl:
            self.ctrlKeyPressed = False
        QtGui.QLineEdit.keyReleaseEvent(self, event)

    def wheelEvent(self, event):
        if self.hasFocus():
            bc = self.base_change_constant
            delta = bc if event.angleDelta().y() > 0 else -bc
            if self.shiftKeyPressed:
                delta = bc*10 if event.angleDelta().y() > 0 else -10*bc
            if self.ctrlKeyPressed:
                delta = bc*0.1 if event.angleDelta().y() > 0 else -0.1*bc
            ctext = self.text()
            try: self.setText('{:.2f}'.format(float(ctext)+delta))
            except: pass
            # fn = self.font()
            # fn.setPointSize(fn.pointSize() +  delta)
            # self.setFont(fn)
            event.accept()
        else: event.ignore()