# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:56:46 2016

@author: sasha
"""

QTVer = 0

def qtVersionToUse():
    try:
        from PyQt4 import QtGui
        QTVer = 4
    except:
        QTVer = 0
        
    if QTVer == 0:
        try:
            from PyQt5 import QtWidgets as QtGui
            QTVer = 5
        except:
            QTVer = 0
            
    return QTVer
    
QTVer = qtVersionToUse()