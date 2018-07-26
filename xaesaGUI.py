# -*- coding: utf-8 -*-
#"""
#Created on Fri Oct  7 12:25:31 2016
#
#@author: sasha
#"""

XAESA_VERSION = "0.01"
GUI_SETTINGS_ID = "XAESA" + XAESA_VERSION

import sys
from sys import exit, argv, version
from os import path, getcwd

import matplotlib.pyplot as plt
from matplotlib import __version__ as mpl_version

__path__=[path.dirname(path.abspath(__file__))]

from .init import QTVer
from .xaesa_exafs_class import xaesa_exafs_class
from .xaesa_xes_class import xaesa_xes_class
from .xaesa_settings import xaesa_settings
from .xaesa_viewer import xaesa_viewer, xaesaViewerWindow

import h5py

def show_exception_and_exit(exc_type, exc_value, tb):
    import traceback
    traceback.print_exception(exc_type, exc_value, tb)
#    raw_input("Press key to exit.")
#    sys.exit(-1)

sys.excepthook = show_exception_and_exit

if QTVer == 4:
    from PyQt4 import QtGui, QtCore
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    print("Python version: " + str(version))
    print("Qt version: " + str(QtCore.qVersion()))
    print("MatPlotLib version: " + str(mpl_version))
    print("h5py version: " + h5py.version.version)
    print("HDF5 version: " + h5py.version.hdf5_version)
    
if QTVer == 5:
    from PyQt5 import QtWidgets as QtGui
    from PyQt5 import QtCore
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    print("Python version: " + str(version))
    print("Qt version: " + str(QtCore.qVersion()))
    print("MatPlotLib version: " + str(mpl_version))
    print("h5py version: " + h5py.version.version)
    print("HDF5 version: " + h5py.version.hdf5_version)


from numpy import asarray, fromstring, genfromtxt, gradient, argmax, zeros, sqrt, sin, transpose, \
                    array, savetxt, copy, delete, arange, exp, argmin, concatenate, all, diff,  \
                    append

#from .ft import FT, BFTWindow, BFT, GETPHASE

#from scipy.interpolate import UnivariateSpline
from scipy.interpolate import Rbf #, UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline

from timeit import default_timer as timer

from .xaesa_deglitch import DGWindow 

from .xaesa_fit import FitWindow
from .xaesa_rdf import RdfWindow, MyStream

from .compare import CompareWindow
from .xaesa_lincombination import LCWindow

#from matplotlib.widgets import RectangleSelector

class MyWindow(QtGui.QMainWindow):
    
#    def resizeEvent(self,resizeEvent):
#
#        print("Information!","Window has been resized...")

    def __init__(self):
        super(MyWindow, self).__init__()
        
        self.exafs_fluo = 0 #0 - exafs, 1 - fluo
        self.calcmju = 1 #calculate mju or take already calculated from file

        self.skiplines = 0
        self.energycol = 0
        self.exafscol = 1
        self.i0col = 7
        self.i1col = 8
        
        self.fluocol = [11,12,13,15,16,17]
        self.setWindowTitle("XAESA - X-ray Absorption and Emission Analytics")
        
        self.dataClasses = []
        
        self.copiedparamsXAS = []
        self.copiedparamsXES = []
        

        self.initUI()
        self.showMaximized()

    def initUI(self):
        
        exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        
        savehdf5Action = QtGui.QAction('Save analysis &to hdf5 ...', self)
        savehdf5Action.setShortcut('Ctrl+S')
        savehdf5Action.triggered.connect(self.savehdf5)
        
        openhdf5Action = QtGui.QAction('Open analysis &from hdf5 ...', self)
        openhdf5Action.setShortcut('Ctrl+A')
        openhdf5Action.triggered.connect(self.openhdf5)
        
        openNexusAction = QtGui.QAction('Open Nexus file ...', self)
        openNexusAction.setShortcut('Ctrl+N')
        openNexusAction.triggered.connect(self.openNexus)
        
        ampPhaSaveAction = QtGui.QAction('Save amplitude and phase ...', self)
        ampPhaSaveAction.triggered.connect(self.ampPhaSave)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(savehdf5Action)
        fileMenu.addAction(openhdf5Action)
        fileMenu.addAction(openNexusAction)
        fileMenu.addAction(exitAction)
        
        toolsMenu = menubar.addMenu('&Tools')
        toolsMenu.addAction(ampPhaSaveAction)
        
        settingsMenu = menubar.addMenu('&Settings')
        self.xaesaSettings = xaesa_settings()
#        self.edtTest = QtGui.QLineEdit('1')
        action = QtGui.QWidgetAction(self)
        action.setDefaultWidget(self.xaesaSettings)
        settingsMenu.addAction(action)
        
        statusBar = self.statusBar()
        self.prgbar = QtGui.QProgressBar()
        self.lblStatus = QtGui.QLabel("Start with opening experimental files or previously saved project.")
        self.lblnpoints = QtGui.QLabel()
        self.lblCurrentFile = QtGui.QLabel()
        statusBar.addWidget(self.lblStatus)
        statusBar.addWidget(self.prgbar)
        statusBar.addWidget(self.lblnpoints)
        statusBar.addWidget(self.lblCurrentFile)

        self.current = -1
        self.currentdir = getcwd()
        self.currenthdf5 = getcwd()
        
        wid = QtGui.QWidget(self)
        self.setCentralWidget(wid)

        #list with spectra
        #from several text files or single hdf5 file
        self.lstSpectra = QtGui.QListWidget()       
        
        self.lstSpectra.itemClicked.connect(self.lstSpectraItemClicked)
        self.lstSpectra.itemDoubleClicked.connect(self.lstSpectraItemDoubleClicked)
        self.lstSpectra.itemActivated.connect(self.lstSpectraItemClicked)
        self.lstSpectra.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
#        self.lstSpectra.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.lstSpectra.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        
        self.lstSpectra.setDefaultDropAction(QtCore.Qt.MoveAction)
        
        self.lstSpectra.installEventFilter(self)
        
#        self.lstSpectra.indexesMoved.connect(self.lstSpectraLayoutChanged)
#        self.lstSpectraModel = self.lstSpectra.model()
#        self.lstSpectraModel.layoutChanged.connect(self.lstSpectraLayoutChanged)
        
        self.lstSpectra.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.lstSpectra.customContextMenuRequested.connect(self.listItemRightClicked)
        
        self.lstSpectraMenu= QtGui.QMenu()

        self.lstSpectraMenu.addAction("Remove Item", self.removefile)        
        self.lstSpectraMenu.addSeparator()        
        self.lstSpectraMenu.addAction("DW factor correction", self.dwcorrection)
        self.lstSpectraMenu.addAction("Shift energy scale by constant", self.shiftenergyscale) 
        self.lstSpectraMenu.addSeparator()        
        self.lstSpectraMenu.addAction("Compare Mju", self.comparemju)        
        self.lstSpectraMenu.addAction("Compare EXAFS", self.compareexafs)        
        self.lstSpectraMenu.addAction("Compare FT", self.compareft)        
        self.lstSpectraMenu.addAction("Compare BFT", self.comparebft)
        self.lstSpectraMenu.addSeparator()  
        self.lstSpectraMenu.addAction("Spectra difference", self.spectradifference)      

        
        #open button
        self.btnOpen = QtGui.QPushButton('Open File(s) ...')
        self.btnRemoveFile = QtGui.QPushButton('Remove File(s) ...')
        self.btnOpen.clicked.connect(self.openfile)
        self.btnRemoveFile.clicked.connect(self.removefile)
        
        self.lblSkipLines = QtGui.QLabel("Skip lines")
        self.lblDataX = QtGui.QLabel("Energy")
        self.lblDataY = QtGui.QLabel("Mju")
        self.lblI0Col = QtGui.QLabel("I0 column")
        self.lblI1Col = QtGui.QLabel("I columns")
#        self.lblFluoCols = QtGui.QLabel("Fluo columns")
        
        self.edtSkipLines = QtGui.QLineEdit(str(self.skiplines))
        self.edtEnergyCol = QtGui.QLineEdit(str(self.energycol))
        self.edtExafsCol = QtGui.QLineEdit(str(self.exafscol))
        self.edtI0Col = QtGui.QLineEdit(str(self.i0col))
        self.edtICol = QtGui.QLineEdit(str(self.i1col))
        self.edtFluoCols = QtGui.QLineEdit("11 12 13 15 16 17")
        self.edtFluoCols.hide()
        
        self.edtSkipLines.setFixedWidth(50)
        self.edtSkipLines.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
        self.edtEnergyCol.setFixedWidth(30)
        self.edtEnergyCol.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
        self.edtExafsCol.setFixedWidth(30)
        self.edtExafsCol.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
        self.edtI0Col.setFixedWidth(30)
        self.edtI0Col.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
#        self.edtI1Col.setFixedWidth(30)
#        self.edtI1Col.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
#        self.edtFluoCols.setFixedWidth(175)
#        self.edtFluoCols.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        
        g1 = QtGui.QHBoxLayout()
        g2 = QtGui.QHBoxLayout()
        g3 = QtGui.QHBoxLayout()
        
        self.gb1 = QtGui.QGroupBox()  
        self.gb2 = QtGui.QGroupBox()
        self.gb3 = QtGui.QGroupBox()
        self.gb1.setLayout(g1)
        self.gb2.setLayout(g2)
        self.gb3.setLayout(g3)
        self.gb1.setFlat(True)
        self.gb2.setFlat(True)
        self.gb3.setFlat(True)
        
        #experiment type
        self.chkDataTypeMju = QtGui.QPushButton("Mju")
        self.chkDataTypeXes = QtGui.QPushButton("XES")
        self.chkDataTypeExafs = QtGui.QPushButton("EXAFS")
        self.chkDataTypeTrans = QtGui.QPushButton("Transient")
        self.chkDataTypeTrans.setVisible(False)
        self.chkDataTypeMju.setCheckable(True)
        self.chkDataTypeXes.setCheckable(True)
        self.chkDataTypeExafs.setCheckable(True)
        self.chkDataTypeTrans.setCheckable(True)
        self.chkDataTypeMju.setAutoExclusive(True)
        self.chkDataTypeXes.setAutoExclusive(True)
        self.chkDataTypeExafs.setAutoExclusive(True)
        self.chkDataTypeTrans.setAutoExclusive(True)
        self.chkDataTypeMju.setStyleSheet("QPushButton:checked { color:green }")
        self.chkDataTypeXes.setStyleSheet("QPushButton:checked { color:green }")
        self.chkDataTypeExafs.setStyleSheet("QPushButton:checked { color:green }")
        self.chkDataTypeTrans.setStyleSheet("QPushButton:checked { color:green }")
        self.chkDataTypeMju.clicked.connect(self.adaptOpeningParams)
        self.chkDataTypeXes.clicked.connect(self.adaptOpeningParams)
        self.chkDataTypeExafs.clicked.connect(self.adaptOpeningParams)
        self.chkDataTypeTrans.clicked.connect(self.adaptOpeningParams)
        

#        self.chkMjuOrExafs = QtGui.QPushButton("Absorption")
#        self.chkMjuOrExafs1 = QtGui.QPushButton("EXAFS")
        self.chkCalcMju = QtGui.QPushButton("Calculate absorption")
        self.chkCalcMju.setCheckable(True)
        self.chkCalcMju.setStyleSheet("QPushButton:checked { color:green }")
        self.chkCalcMju.clicked.connect(self.adaptOpeningParams)
        self.chkTransOrFluo = QtGui.QPushButton("Transmission")
        self.chkTransOrFluo1 = QtGui.QPushButton("Fluorescence")
        self.chkTransOrFluo.clicked.connect(self.adaptOpeningParams)
        self.chkTransOrFluo1.clicked.connect(self.adaptOpeningParams)

        self.chkTransOrFluo.setCheckable(True)
        self.chkTransOrFluo1.setCheckable(True)
#        self.chkMjuOrExafs.setCheckable(True)
#        self.chkMjuOrExafs1.setCheckable(True)
        
        self.chkTransOrFluo.setAutoExclusive(True)
        self.chkTransOrFluo1.setAutoExclusive(True)
#        self.chkMjuOrExafs.setAutoExclusive(True)
#        self.chkMjuOrExafs1.setAutoExclusive(True)
        
        self.chkTransOrFluo.setStyleSheet("QPushButton:checked { color:green }")
        self.chkTransOrFluo1.setStyleSheet("QPushButton:checked { color:green }")
#        self.chkMjuOrExafs.setStyleSheet("QPushButton:checked { color:green }")
#        self.chkMjuOrExafs1.setStyleSheet("QPushButton:checked { color:green }")
        
        
        self.chkTransOrFluo.setChecked(True)
        self.chkDataTypeMju.setChecked(True) 

        g1.addWidget(self.chkDataTypeMju)
        g1.addWidget(self.chkDataTypeXes)
        g1.addWidget(self.chkDataTypeExafs)
        g1.addWidget(self.chkDataTypeTrans)
#        g1.addWidget(self.lblEnergyCol)
#        g1.addWidget(self.edtEnergyCol)
#        g1.addWidget(self.lblExafsCol)
#        g1.addWidget(self.edtExafsCol) 
        
        g2.addWidget(self.chkCalcMju)
#        g2.addWidget(self.lblI0Col)
#        g2.addWidget(self.edtI0Col)
#        g2.addWidget(self.lblI1Col)
#        g2.addWidget(self.edtI1Col)

        g3.addWidget(self.chkTransOrFluo)
        g3.addWidget(self.chkTransOrFluo1)
#        g3.addWidget(self.lblFluoCols)
#        g3.addWidget(self.edtFluoCols)

        lopenparams = QtGui.QGridLayout()
        lopenparams.addWidget(self.lblSkipLines, 0, 0)
        lopenparams.addWidget(self.edtSkipLines, 0, 1)
        
        lopenparams.addWidget(self.gb1, 1, 0, 1, 6)
        
        lopenparams.addWidget(self.lblDataX, 2, 0)
        lopenparams.addWidget(self.edtEnergyCol, 2, 1)
        lopenparams.addWidget(self.lblDataY, 2, 2)
        lopenparams.addWidget(self.edtExafsCol, 2, 3)

        lopenparams.addWidget(self.gb2, 3, 0, 1, 2)
        lopenparams.addWidget(self.gb3, 3, 3, 1, 3)
        
        lopenparams.addWidget(self.lblI0Col, 4, 0)
        lopenparams.addWidget(self.edtI0Col, 4, 1)
        lopenparams.addWidget(self.lblI1Col, 4, 2)        
        lopenparams.addWidget(self.edtICol, 4, 3)
        lopenparams.addWidget(self.edtFluoCols, 4, 3)
        
        
#        lopenparams.addWidget(self.lblFluoCols, 3,2)
#        lopenparams.addWidget(self.edtFluoCols, 3,3, 1, 3)

        lh = QtGui.QHBoxLayout()
        
        l1 = QtGui.QVBoxLayout()
        lh.addWidget(self.btnOpen)
        lh.addWidget(self.btnRemoveFile)
        l1.addLayout(lh)
        l1.addLayout(lopenparams)
        l1.addWidget(self.lstSpectra)

    #####################         XAS FIGURES
        #Figures 
        self.fig = plt.figure(0)
        self.ax_abs = plt.subplot2grid((2,3), (0,0))
        
        #init lines
        self.abs_scatter, = self.ax_abs.plot([], [], 'o', picker=True)
        self.abs_scatter.set_markersize(2)

        self.mjub_line, self.mju0_mjub_line = self.ax_abs.plot([], [], [], [])

        self.mjub_line.set_color('g')
        self.mju0_mjub_line.set_color('g')
        
        self.ax_abs2 = self.ax_abs.twinx()
        
        self.deriv_line, = self.ax_abs2.plot([],[])
        self.deriv_line.set_color('b')
        self.deriv_line.set_alpha(0.3)
        
        self.lineE0 = self.ax_abs.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE1 = self.ax_abs.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE2 = self.ax_abs.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE3 = self.ax_abs.axvline(0, color='g', linestyle='--', lw=1)
        
        self.ax_abs.set_xlabel('Energy, eV')
        self.ax_abs.set_ylabel('Absorption, a.u.')
        
        
        self.ax_exafs = plt.subplot2grid((2,3), (1,0))
        

        self.exafs_line, = self.ax_exafs.plot([], [], label = "EXAFS")
        self.exafs_line.set_color('k')  
        self.exafssm_line, = self.ax_exafs.plot([], [], label = "zlc")
        self.exafssm_line.set_color('r')  
        
        self.ax_exafs_legend = self.ax_exafs.legend(handles=[self.exafs_line, self.exafssm_line])
        
        self.ax_exafs.axhline(y=0, linewidth=1, color = 'k')

        self.ax_exafssm = plt.subplot2grid((2,3), (0,1), colspan=2)
        
        self.ax_exafssm.axhline(y=0, linewidth=1, color = 'k')
        
        self.exafsd_line, = self.ax_exafssm.plot([], [], label = "EXAFS zlc")
        self.bftexafs_line, = self.ax_exafssm.plot( [], [], label = "BFT")
        self.bftexafs_line.set_color('r')
        
        
        self.ax_exafs.set_xlabel('Wavevector k, $\AA^{-1}$')
        self.ax_exafs.set_ylabel('EXAFS, $\AA^{-2}$')
        self.ax_exafssm.set_xlabel('Wavevector k, $\AA^{-1}$')
        self.ax_exafssm.set_ylabel('EXAFS, $\AA^{-2}$')
        
        self.ax_exafsm_legend = self.ax_exafssm.legend(handles=[self.exafsd_line, self.bftexafs_line])
        
        
        self.ax_ft = plt.subplot2grid((2,3), (1,1), colspan=2)
        
        
        self.efr_line, self.efi_line = self.ax_ft.plot( [],  [], [],  [], label = "FT")
        self.efr_line.set_color('g')
        self.efi_line.set_color('g')
        self.efi_line.set_linestyle('dotted')
        self.efrsm_line, self.efism_line = self.ax_ft.plot([],  [], [],  [], label = "FT zlc")
        self.efrsm_line.set_color('b')
        self.efism_line.set_color('b')
        self.efism_line.set_linestyle('dotted')
        self.bftwinr_line, self.bftwini_line = self.ax_ft.plot([], [], [], [], label = "BFT FT")
        self.bftwinr_line.set_color('r')
        self.bftwini_line.set_color('r')
        self.bftwini_line.set_linestyle('dotted')
        
        self.ax_ft_legend = self.ax_ft.legend(handles=[self.efr_line, self.efrsm_line, self.bftwinr_line])
        
        self.canv = FigureCanvas(self.fig)
        self.tbar = NavigationToolbar(self.canv, self)
        fnt = self.tbar.font()
        fnt.setPointSize(20)
        self.tbar.setFont(fnt)
        
        self.ax_ft.set_xlabel('Distance R, $\AA$')
        self.ax_ft.set_ylabel('Fourier transform, $\AA^{-3}$')
        
        self.fig.tight_layout()
        
        self.fig.canvas.mpl_connect('pick_event', self.onpick3)
        
        self.ax_ft.axhline(y=0, linewidth=1, color = 'k')

        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbar)
        lfig.addWidget(self.canv)
        
        #Frame for xas figures
        self.frameXasFig = QtGui.QFrame()
        self.frameXasFig.setFrameShape(QtGui.QFrame.Panel)
        self.frameXasFig.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameXasFig.setLayout(lfig)
        self.frameXasFig.hide()
        
#XAS FIGURES END
        
    ##################             XES FIGURES

        self.figXes = plt.figure(11)
        self.ax_xes = plt.subplot2grid((2,2), (0,0))
        
        #init lines
        #self.abs_scatter = self.ax_abs.scatter([], [], s=1, picker=True)
        self.lineXesOriginal, = self.ax_xes.plot([], [], 'o', picker=True)
        self.lineXesOriginal.set_markersize(2)

        self.lineXesBkgr, = self.ax_xes.plot([], [])
        self.lineXesBkgr.set_color('g')
        
        self.lineE0xes = self.ax_xes.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE1xes = self.ax_xes.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE2xes = self.ax_xes.axvline(0, color='g', linestyle='--', lw=1)
        self.lineE3xes = self.ax_xes.axvline(0, color='g', linestyle='--', lw=1)
        
        self.lineEminAnorm = self.ax_xes.axvline(0, color='r', linestyle='--', lw=1)
        self.lineEmaxAnorm = self.ax_xes.axvline(0, color='r', linestyle='--', lw=1)

        
        self.ax_xes.set_xlabel('Energy, eV')
        self.ax_xes.set_ylabel('Emission intensity, a.u.')
        
        
        self.ax_xes_bc = plt.subplot2grid((2,2), (0,1))
        

        self.lineXesBkgrCorr, = self.ax_xes_bc.plot([], []) #, label = "EXAFS")
        self.lineXesBkgrCorr.set_color('k')  
#        self.exafssm_line, = self.ax_xes_bc.plot([], []) #, label = "zlc")
#        self.exafssm_line.set_color('r')  
        
#        self.ax_exafs.legend(handles=[self.exafs_line, self.exafssm_line])
        
        self.ax_xes_bc.axhline(y=0, linewidth=1, color = 'k')

        self.ax_xes_norm = plt.subplot2grid((2,2), (1,0), colspan=2)
        
        self.ax_xes_norm.axhline(y=0, linewidth=1, color = 'k')
        
        self.lineXesNorm, = self.ax_xes_norm.plot([], []) #, label = "EXAFS zlc")
        
        self.canvXes = FigureCanvas(self.figXes)
        self.tbarXes = NavigationToolbar(self.canvXes, self)
        
#        self.ax_ft.set_xlabel('Distance R, $\AA$')
#        self.ax_ft.set_ylabel('Fourier transform, $\AA^{-3}$')
        
        self.figXes.tight_layout()
        
#        self.figXes.canvas.mpl_connect('pick_event', self.onpick3)
        
#        self.ax_ft.axhline(y=0, linewidth=1, color = 'k')
    
#        
    

        lfig = QtGui.QVBoxLayout()
        lfig.addWidget(self.tbarXes)
        lfig.addWidget(self.canvXes)
        
        #Frame for XES figures
        self.frameXesFig = QtGui.QFrame()
        self.frameXesFig.setFrameShape(QtGui.QFrame.Panel)
        self.frameXesFig.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameXesFig.setLayout(lfig)
#        self.frameXesFig.hide()
        
#XES FIGURES END
        
    ######################   XES parameters
        lblE0xes = QtGui.QLabel("E0")
        lblE1xes = QtGui.QLabel("E1")
        lblE2xes = QtGui.QLabel("E2")
        lblE3xes = QtGui.QLabel("E3")
        lblEMinANormxes = QtGui.QLabel("Emin for area norm")
        lblEAMaxNormxes = QtGui.QLabel("Emax for area norm")
        
        self.edtXesE0 = QtGui.QLineEdit()
        self.edtXesE1 = QtGui.QLineEdit()
        self.edtXesE2 = QtGui.QLineEdit()
        self.edtXesE3 = QtGui.QLineEdit()
        
        self.edtXesEANormMin = QtGui.QLineEdit()
        self.edtXesEANormMax = QtGui.QLineEdit()
        
        self.edtXesE0.returnPressed.connect(self.xesExtractRedo)
        self.edtXesE1.returnPressed.connect(self.xesExtractRedo)
        self.edtXesE2.returnPressed.connect(self.xesExtractRedo)
        self.edtXesE3.returnPressed.connect(self.xesExtractRedo)
        
        self.edtXesEANormMin.returnPressed.connect(self.xesExtractRedo)
        self.edtXesEANormMax.returnPressed.connect(self.xesExtractRedo)
        
        lparams = QtGui.QGridLayout()
        lparams.addWidget(lblE0xes, 0, 0)
        lparams.addWidget(lblE1xes, 1, 0)
        lparams.addWidget(lblE2xes, 2, 0)
        lparams.addWidget(lblE3xes, 3, 0)
        
        lparams.addWidget(self.edtXesE0, 0, 1)
        lparams.addWidget(self.edtXesE1, 1, 1)
        lparams.addWidget(self.edtXesE2, 2, 1)
        lparams.addWidget(self.edtXesE3, 3, 1)
        
        lparams.addWidget(lblEMinANormxes, 0, 2)
        lparams.addWidget(lblEAMaxNormxes, 1, 2)
        
        lparams.addWidget(self.edtXesEANormMin, 0, 3)
        lparams.addWidget(self.edtXesEANormMax, 1, 3)
        
        
        self.frameXesParams = QtGui.QFrame()
        self.frameXesParams.setFrameShape(QtGui.QFrame.Panel)
        self.frameXesParams.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameXesParams.setLayout(lparams)
        self.frameXesParams.hide()
        
        
#XES parameters END     

        #Energy params
        self.lblE0 = QtGui.QLabel("E0")
        self.lblE1 = QtGui.QLabel("E1")
        self.lblE2 = QtGui.QLabel("E2")
        self.lblE3 = QtGui.QLabel("E3")
        self.chkExafsNormalization = QtGui.QCheckBox("Normalize by edge at energy")
        self.lblkpow = QtGui.QLabel("k power")
        self.lblsm = QtGui.QLabel("3rd zero-line correction (zlc)")
        self.lblmju0poldegree = QtGui.QLabel("Mju0 polynomial degree")
        self.edtE0 = QtGui.QLineEdit()
        self.edtE1 = QtGui.QLineEdit()
        self.edtE2 = QtGui.QLineEdit()
        self.edtE3 = QtGui.QLineEdit()
        self.edtExafsNormEnergy = QtGui.QLineEdit()
        self.edtkpow = QtGui.QLineEdit("2")
        self.edtsm = QtGui.QLineEdit("0")
        self.edtmju0poldegree = QtGui.QLineEdit("4")
        self.btnDeglitching = QtGui.QPushButton('Remove glitches ...')
        self.btnDeglitching.clicked.connect(self.removeglitches)
        
        
        lparams = QtGui.QGridLayout()
        lparams1 = QtGui.QHBoxLayout()
        lparams1.addWidget(self.lblE0)
        lparams1.addWidget(self.edtE0)
        lparams1.addWidget(self.lblE1)
        lparams1.addWidget(self.edtE1)
        lparams1.addWidget(self.lblE2)
        lparams1.addWidget(self.edtE2)
        lparams1.addWidget(self.lblE3)
        lparams1.addWidget(self.edtE3)
        lparams.addLayout(lparams1, 0,0,1,4)
        lparams.addWidget(self.chkExafsNormalization, 1, 0)
        lparams.addWidget(self.edtExafsNormEnergy, 1, 1)
#        lparams.addWidget(self.lblE0, 0, 0)
#        lparams.addWidget(self.lblE1, 0, 1)
#        lparams.addWidget(self.lblE2, 0, 2)
#        lparams.addWidget(self.lblE3, 0, 3)
        lparams.addWidget(self.lblkpow, 2, 2)
        
        
        
        
#        lparams.addWidget(self.edtE0, 1, 0)
#        lparams.addWidget(self.edtE1, 1, 1)
#        lparams.addWidget(self.edtE2, 1, 2)
#        lparams.addWidget(self.edtE3, 1, 3)
        lparams.addWidget(self.edtkpow, 2, 3)
        lparams.addWidget(self.lblsm, 2, 0)
        lparams.addWidget(self.edtsm, 2, 1)
        lparams.addWidget(self.lblmju0poldegree, 3, 0)
        lparams.addWidget(self.edtmju0poldegree, 3, 1)
        lparams.addWidget(self.btnDeglitching, 3, 2, 1, 2)
        self.edtE0.returnPressed.connect(self.extractParamsChange)
        self.edtE1.returnPressed.connect(self.extractParamsChange)
        self.edtE2.returnPressed.connect(self.extractParamsChange)
        self.edtE3.returnPressed.connect(self.extractParamsChange)
        self.edtkpow.returnPressed.connect(self.kPowerChange)
        self.edtsm.returnPressed.connect(self.extractParamsChange)
        self.edtmju0poldegree.returnPressed.connect(self.extractParamsChange)
        self.edtExafsNormEnergy.returnPressed.connect(self.extractParamsChange)
        
        #Frame for lparams
        self.frameLParams = QtGui.QFrame()
        self.frameLParams.setFrameShape(QtGui.QFrame.Panel)
        self.frameLParams.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameLParams.setLayout(lparams)

        
        
        self.lblkmin = QtGui.QLabel("K min")
        self.lblkmax = QtGui.QLabel("K max")
        self.lbldk = QtGui.QLabel("dK")
        self.lblrmin = QtGui.QLabel("R min")
        self.lblrmax = QtGui.QLabel("R max")
        self.lbldr = QtGui.QLabel("dR")

        self.edtkmin = QtGui.QLineEdit("0.5")
        self.edtkmax = QtGui.QLineEdit("18")
        self.edtdk = QtGui.QLineEdit("0.05")
        self.edtrmin = QtGui.QLineEdit("0")
        self.edtrmax = QtGui.QLineEdit("6")
        self.edtdr = QtGui.QLineEdit("0.02")
        self.edtkmin.returnPressed.connect(self.ftBftParamsChange)
        self.edtkmax.returnPressed.connect(self.ftBftParamsChange)
        self.edtdk.returnPressed.connect(self.ftBftParamsChange)
        self.edtrmin.returnPressed.connect(self.ftBftParamsChange)
        self.edtrmax.returnPressed.connect(self.ftBftParamsChange)
        self.edtdr.returnPressed.connect(self.ftBftParamsChange)

        lftparams = QtGui.QGridLayout()
        lftparams.addWidget(self.lblkmin, 0, 0)
        lftparams.addWidget(self.lblkmax, 0, 1)
        lftparams.addWidget(self.lbldk,   0, 2)
        lftparams.addWidget(self.edtkmin, 1, 0)
        lftparams.addWidget(self.edtkmax, 1, 1)
        lftparams.addWidget(self.edtdk,   1, 2)

        lftparams.addWidget(self.lblrmin, 2, 0)
        lftparams.addWidget(self.lblrmax, 2, 1)
        lftparams.addWidget(self.lbldr,   2, 2)
        lftparams.addWidget(self.edtrmin, 3, 0)
        lftparams.addWidget(self.edtrmax, 3, 1)
        lftparams.addWidget(self.edtdr,   3, 2)
        
        #Frame for lftparams
        self.frameLFTParams = QtGui.QFrame()
        self.frameLFTParams.setFrameShape(QtGui.QFrame.Panel)
        self.frameLFTParams.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameLFTParams.setLayout(lftparams)

        self.lblrminbft = QtGui.QLabel("BFT R min")
        self.lblrmaxbft = QtGui.QLabel("BFT R max")
        self.lblbftwindowparam = QtGui.QLabel("BFT window parameter")

        self.edtrminbft = QtGui.QLineEdit("0")
        self.edtrmaxbft = QtGui.QLineEdit("6")
        self.edtbftwindowparam = QtGui.QLineEdit("0.1")
        self.edtrminbft.returnPressed.connect(self.bftParamsChange)
        self.edtrmaxbft.returnPressed.connect(self.bftParamsChange)
        self.edtbftwindowparam.returnPressed.connect(self.bftParamsChange)
        self.btnFit = QtGui.QPushButton('Fit 1st sphere ...')
        self.btnFit.clicked.connect(self.fit)
        self.btnRdf = QtGui.QPushButton('Reconstruct RDF ...')
        self.btnRdf.clicked.connect(self.rdf)


        lbftparams = QtGui.QGridLayout()
        lbftparams.addWidget(self.lblrminbft, 0,0)
        lbftparams.addWidget(self.lblrmaxbft, 0,1)
        lbftparams.addWidget(self.lblbftwindowparam, 0,2)

        lbftparams.addWidget(self.edtrminbft, 1,0)
        lbftparams.addWidget(self.edtrmaxbft, 1,1)
        lbftparams.addWidget(self.edtbftwindowparam, 1,2)
        lbftparams.addWidget(self.btnFit, 2,0, 1, 1)
        lbftparams.addWidget(self.btnRdf, 2,1, 1, 2)
        
        #Frame for lbftparams
        self.frameLBFTParams = QtGui.QFrame()
        self.frameLBFTParams.setFrameShape(QtGui.QFrame.Panel)
        self.frameLBFTParams.setFrameShadow(QtGui.QFrame.Sunken)
        self.frameLBFTParams.setLayout(lbftparams)
        
        #compare buttons
        self.mnuCompareXAS = QtGui.QMenu()
        self.mnuCompareXAS.addAction('Compare mju...', self.comparemju)
        self.mnuCompareXAS.addAction('Compare XANES...', self.comparexanes)
        self.mnuCompareXAS.addAction('Compare EXAFS...', self.compareexafs)
        self.mnuCompareXAS.addAction('Compare FT...', self.compareft)
        self.mnuCompareXAS.addAction('Compare BFT...', self.comparebft)
        self.mnuCompareXAS.addSeparator()
        self.mnuCompareXAS.addAction('Compare i0', lambda: self.compare('i0'))
        self.mnuCompareXAS.addAction('Compare i1', lambda: self.compare('i1'))
        self.mnuCompareXAS.addAction('Compare i fluorescence', lambda: self.compare('ifluo'))
        self.mnuCompareXAS.addSeparator()
        self.mnuCompareXAS.addAction('Compare amplitude', lambda: self.compare('amp'))
        self.mnuCompareXAS.addAction('Compare phase', lambda: self.compare('pha'))
        self.mnuCompareXAS.addSeparator()
        self.mnuCompareXAS.addAction('Compare mju rebin vs original', lambda: self.compare('mjuRebinVSoriginal'))
#        
        self.mnuCompareXES = QtGui.QMenu()
        self.mnuCompareXES.addAction('Compare XES...', self.compareXes)
        self.mnuCompareXES.addAction('Compare XES area normalized...', self.compareXesAnorm)

        
        self.btnCompareXas = QtGui.QPushButton('Compare ...')
#        self.btnCompareXas.setMenu(self.mnuCompareXAS)
        
        #save buttons
        self.mnuSaveXAS = QtGui.QMenu()
        self.mnuSaveXAS.addAction('Save mju...', self.savemju)
        self.mnuSaveXAS.addAction('Save XANES...', self.savexanes)
        self.mnuSaveXAS.addAction('Save EXAFS...', self.saveexafs)
        self.mnuSaveXAS.addAction('Save FT...', self.saveft)
        self.mnuSaveXAS.addAction('Save BFT...', self.savebft)
        
        self.mnuSaveXES = QtGui.QMenu()
        self.mnuSaveXES.addAction('Save XES...', self.saveXes)
        self.mnuSaveXES.addAction('Save XES area normalized...', self.saveXesAnorm)

        self.btnSaveXas = QtGui.QPushButton('Save ...')
#        self.btnSaveXas.setMenu(self.mnuSaveXAS)
        
        lh = QtGui.QHBoxLayout()
        lh.addWidget(self.btnCompareXas)
        lh.addWidget(self.btnSaveXas)
        
        l1.addLayout(lh)

        self.chkSaveinOneFile = QtGui.QCheckBox("Save selected in single file")        
        self.chkSaveinOneFile.setChecked(True)
    
        
        l1.addLayout(lh)
        l1.addWidget(self.chkSaveinOneFile)


        self.btnAverage = QtGui.QPushButton('Average')
        
        self.mnuAverageXAS = QtGui.QMenu()
        self.mnuAverageXAS.addAction('Average mju', self.averageMju)
        self.mnuAverageXAS.addAction('Average EXAFS', self.averageExafs)
        
        self.mnuAverageXES = QtGui.QMenu()
        self.mnuAverageXES.addAction('Average XES', self.averageXes)
        

        self.btnRebin = QtGui.QPushButton('Rebin')
        
        self.edtRbfSmooth = QtGui.QLineEdit('1')
        action = QtGui.QWidgetAction(self)
        action.setDefaultWidget(self.edtRbfSmooth)
        
        self.mnuRebinXAS = QtGui.QMenu()
        self.mnuRebinXAS.addAction('Rebin mju', self.rebinMju)
        self.mnuRebinXAS.addSeparator()
        self.mnuRebinXAS.addAction(action)
        self.mnuRebinXAS.addAction('Rebin mju-RBF-smooth', self.rebinMjuRbfSmooth)
        self.mnuRebinXAS.addSeparator()
        self.mnuRebinXAS.addAction('Rebin EXAFS', self.rebinExafs)
        self.mnuRebinXAS.addAction('Return to original mju', self.backToOriginalMju)
        

        
        
        self.mnuRebinXES = QtGui.QMenu()
        self.mnuRebinXES.addAction('Rebin XES', self.rebinXes)

        self.btnLinCombination = QtGui.QPushButton('Linear combination...')
        self.btnLinCombination.clicked.connect(self.LinCombination)
#        self.chkPicker = QtGui.QCheckBox("Allow to remove experimental points")
#        self.chkPicker.stateChanged.connect(self.redraw_for_remove)

        self.btnCopyParams = QtGui.QPushButton('Copy params')
        self.btnCopyParams.clicked.connect(self.CopyParams)
        self.btnApplytoSelected = QtGui.QPushButton('Apply params to selected')
        self.btnApplytoSelected.clicked.connect(self.ApplytoSelected)
        lh1 = QtGui.QHBoxLayout()
        lh1.addWidget(self.btnCopyParams)
        lh1.addWidget(self.btnApplytoSelected)

        lbtn = QtGui.QGridLayout()
        
        lbtn.addLayout(lh1, 0, 0, 1, 2)
        
        lbtn.addWidget(self.btnAverage, 1, 0)
        lbtn.addWidget(self.btnRebin, 1, 1)
        lbtn.addWidget(self.btnLinCombination, 2, 0, 1, 2)

        layout = QtGui.QGridLayout()
        layout.addLayout(l1,      0, 0)
        layout.addWidget(self.frameXasFig,    0, 1, 1,3)
        layout.addWidget(self.frameXesFig,    0, 1, 1,3)
        

        layout.addLayout(lbtn, 1, 0)
        layout.addWidget(self.frameLParams, 1, 1)
        layout.addWidget(self.frameXesParams, 1, 1)        
        layout.addWidget(self.frameLFTParams, 1, 2)
        layout.addWidget(self.frameLBFTParams, 1, 3)
        
        layout.setColumnStretch(0,1)
        layout.setColumnStretch(1,1)
        layout.setColumnStretch(2,1)
        layout.setColumnStretch(3,1)
        
        lftparams.setAlignment(QtCore.Qt.AlignTop)  
        lbftparams.setAlignment(QtCore.Qt.AlignTop)

        wid.setLayout(layout)
        
        self.initialdir = ""
        
        self.fig.tight_layout()
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.canv.draw()
        
        self.adaptOpeningParams()
        
        self.leg_lines = self.ax_exafsm_legend.get_lines() + \
                    self.ax_ft_legend.get_lines() + self.ax_ft_legend.get_lines()
        print(self.leg_lines)
#        leg_lines.append(self.ax_exafs_legend.get_lines())
#        leg_lines.append(self.ax_exafsm_legend.get_lines())
#        leg_lines.append(self.ax_ft_legend.get_lines())
        
        self.lines = [self.exafsd_line, self.bftexafs_line,
                 self.efr_line,  self.efrsm_line, self.bftwinr_line,
                 self.efi_line,  self.efism_line,  self.bftwini_line]
        
        print(self.lines)
        
        
        # we will set up a dict mapping legend line to orig line, and enable
        # picking on the legend line
        self.linedXas = dict()
        for leglineXas, origlineXas in zip(self.leg_lines, self.lines):
            leglineXas.set_picker(5)  # 5 pts tolerance
            self.linedXas[leglineXas] = origlineXas

#        if self.mode == 1:
#            self.lined1 = dict()
#            for legline, origline in zip(leg_lines, self.lines1):
##                legline.set_picker(5)  # 5 pts tolerance
#                self.lined1[legline] = origline

        self.show()

        
    def openaddfiledialog_qtgui(self):
        
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["XDI files (*.xdi)",
                            "FIO files (*.fio)", 
                            "Text files (*.txt)", 
                            "DAT files (*.dat)", 
                            "All files (*.*)"])
#        dlg.setDirectory(self.currentdir)
#        filenames = QStringList()
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames
        else:
            return []
        
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
    
    def directoryfiledialog_qtgui(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.Directory)
#        dlg.setDirectory(self.currentdir)
        dlg.setOption(QtGui.QFileDialog.ShowDirsOnly, True)
        if dlg.exec_():
            dlist =  dlg.selectedFiles()
            return dlist[0]
        else:
            return ""

    def openfile(self):
        
#        #apply file format settings
        self.skiplines = int(self.edtSkipLines.text())
        self.i0col = int(self.edtI0Col.text())
        self.i1col = int(self.edtICol.text())
        self.energycol = int(self.edtEnergyCol.text())
        self.exafscol = int(self.edtExafsCol.text())
        self.fluocol = fromstring(self.edtFluoCols.text(), dtype=int, sep=' ')
        
        self.exafs_fluo = int(not self.chkTransOrFluo.isChecked())
        self.calcmju = int(self.chkCalcMju.isChecked())
#        

        self.fn = self.openaddfiledialog_qtgui()
        if self.fn == []:
            self.lblStatus.setText("No files opened")
            return
        
        self.fn.sort()
        
        fn1 = []
        
        for i in range(0, len(self.fn)):
            head, tail = path.split(str(self.fn[i]))
            fn1.append(tail)
            
        self.currentdir = head        
        
        if self.chkDataTypeExafs.isChecked(): #only exafs
            msgBox = QtGui.QMessageBox()
            msgBox.setText('Select k power used?')
            msgBox.addButton(QtGui.QPushButton('0'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('1'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('2'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('3'), QtGui.QMessageBox.AcceptRole)
            kpower = msgBox.exec_()
            
        startTime = timer()

        for i in range(0, len(self.fn)):
            
            #test code for automatic determination of the file structure
            try:
                hl, fl, n_hashtag, ncol, npoints = self.openaddfile_asciistructure(self.fn[i])
                #check open settings
                
            except:
                pass
            
            self.lstSpectra.addItem(fn1[i])
    ############### XES   ###################################################################################
            if self.chkDataTypeXes.isChecked(): #XES data
                self.dataClasses.append(xaesa_xes_class())
                self.dataClasses[-1].name = fn1[i]
                try:
                    self.dataClasses[-1].energy, self.dataClasses[-1].xes = genfromtxt(str(self.fn[i]),
                                                        comments = "#",
                                                        skip_header = self.skiplines, 
                                                        usecols=(self.energycol, self.exafscol), 
                                                        unpack=True)
                except:
                    self.openaddexception_cleaner(self.fn[i])
                    return

                self.dataClasses[-1].E0 = self.dataClasses[-1].energy[0]
                self.dataClasses[-1].E1 = self.dataClasses[-1].energy[int(len(self.dataClasses[-1].energy)/5)]
                self.dataClasses[-1].E2 = self.dataClasses[-1].energy[int(4*len(self.dataClasses[-1].energy)/5)]
                self.dataClasses[-1].E3 = self.dataClasses[-1].energy[len(self.dataClasses[-1].energy)-1]
                self.dataClasses[-1].eAreaNormMin = self.dataClasses[-1].energy[0]
                self.dataClasses[-1].eAreaNormMax = self.dataClasses[-1].energy[len(self.dataClasses[-1].energy)-1]
                
                #process opened file
                self.dataClasses[-1].removeBackground()
                self.dataClasses[-1].areaNormalize()
                
                sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                sel_item.setForeground(QtCore.Qt.darkMagenta)
            
            
    ############### EXAFS   ###################################################################################
            if self.chkDataTypeExafs.isChecked(): #only exafs
                self.dataClasses.append(xaesa_exafs_class(3))
                self.dataClasses[-1].name = fn1[i]
                try:
                    self.dataClasses[-1].k, self.dataClasses[-1].exafs = genfromtxt(str(self.fn[i]),
                                                        comments = "#",
                                                        skip_header = self.skiplines, 
                                                        usecols=(self.energycol, self.exafscol), 
                                                        unpack=True)
                except:
                    self.openaddexception_cleaner(self.fn[i])
                    return
                
                    
                self.dataClasses[-1].E0 = 0
                self.dataClasses[-1].E1 = 0
                self.dataClasses[-1].E2 = 0
                self.dataClasses[-1].E3 = 0
    
                self.dataClasses[-1].kPower = kpower
                
                self.dataClasses[-1].processExpData()
                
                sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                sel_item.setForeground(QtCore.Qt.blue)

    ############### Calculate Mju   #########
            if self.chkDataTypeMju.isChecked(): #mju
                if self.chkCalcMju.isChecked():
                    if self.exafs_fluo == 0: # xas
                        self.dataClasses.append(xaesa_exafs_class(0))
                        self.dataClasses[-1].name = fn1[i]
                        try:
                            self.dataClasses[-1].energy, \
                            self.dataClasses[-1].i0, \
                            self.dataClasses[-1].i1 = genfromtxt(str(self.fn[i]),
                                                        comments = "#",
                                                        skip_header = self.skiplines, 
                                                        usecols=(self.energycol, self.i0col, self.i1col), 
                                                        unpack=True)
                        except:
                            print("error opening file")

                        try:
                            self.dataClasses[-1].processExpData()
                        except:
                            print("error processing data")
                            continue
                        sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                        sel_item.setForeground(QtCore.Qt.darkGreen)
                        
#                        if len(self.dataClasses[-1].energy) > 2000: #make auto rebin
#                            self.dataClasses[-1].changeToRebinedMju()
#                            sel_item.setText( sel_item.text() + "(rebin)")
                        
                    if self.exafs_fluo == 1: #fluo
                        self.dataClasses.append(xaesa_exafs_class(1))
                        self.dataClasses[-1].name = fn1[i]
                        
                        self.dataClasses[-1].energy, \
                        self.dataClasses[-1].i0 = genfromtxt(str(self.fn[i]), comments = "#", skip_header = self.skiplines, usecols=(self.energycol, self.i0col), unpack=True)
                        
                        self.dataClasses[-1].ifluo = genfromtxt(str(self.fn[i]), comments = "#", skip_header = self.skiplines, usecols=self.fluocol, unpack=True)
                        
                        if self.dataClasses[-1].ifluo.ndim == 1:
                            self.dataClasses[-1].ifluo = asarray([self.dataClasses[-1].ifluo])
                            
                        print(self.dataClasses[-1].ifluo)
                        try:    
                            self.dataClasses[-1].processExpData()
                        except:
                            print("error processing data")
                            continue
                        
                        sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                        sel_item.setForeground(QtCore.Qt.darkGreen)
                        
#                        if len(self.dataClasses[-1].energy) > 2000: #make auto rebin
#                            self.dataClasses[-1].changeToRebinedMju()
#                            sel_item.setText( sel_item.text() + "(rebin)")
    ###############  Mju from file  ########               
                else:
                    self.dataClasses.append(xaesa_exafs_class(2)) # mju
                    self.dataClasses[-1].name = fn1[i]
                    try:
                        self.dataClasses[-1].energy, \
                        self.dataClasses[-1].mju = genfromtxt(  str(self.fn[i]), comments = "#", 
                                                                skip_header = self.skiplines, 
                                                                usecols=(self.energycol, self.exafscol), 
                                                                unpack=True)
                    except:
                        print("error opening mju")
                        pass
                    try:   
                        self.dataClasses[-1].processExpData()
                    except:
                        print("error processing data")
                        continue
                    sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                    sel_item.setForeground(QtCore.Qt.darkGreen)
                    
#                    if len(self.dataClasses[-1].energy) > 2000: #make auto rebin
#                        self.dataClasses[-1].changeToRebinedMju()
#                        sel_item.setText( sel_item.text() + "(rebin)")
        
        endTime = timer()
        self.lblStatus.setText("{:d} files opened in {:.4f} seconds".format(len(self.fn), endTime - startTime) )
#        self.lblStatus.setText("Files opened.")    
        
    def openaddfile_asciistructure(self, filename):
    #test code for automatic determination of the file structure
        ncol = []
        f = open(filename)
#        data = f.read()
        nlines = 0
        n_hashtag = 0
        for line in f:
            cols = line.split()
            ncol.append( len(cols))
            nlines = nlines+1
            if line[0] == '#':
                n_hashtag = n_hashtag + 1
        curr = ncol[0]
        starti = [0]
        leng = [1]
        for j in range(1, len(ncol)):
            if curr==ncol[j]:
                leng[-1] = leng[-1]+1
            else:
                curr = ncol[j]
                starti.append(j)
                leng.append(1)
        hl = starti[ argmax(leng) ]
        fl = nlines - hl - leng[ argmax(leng) ]
        print("File ", filename)
        print("Header lines : ", hl)
        print("Footer lines : ", fl)
        print("Lines started with # : ", n_hashtag)
        print("Number of columns : ", ncol[hl+1])
        print("Number of data points : ", leng[ argmax(leng)] )
        print("")
#        print(starti)
#        print(leng)
        f.close()
        return hl, fl, n_hashtag, ncol[hl+1], leng[ argmax(leng)]
        
    def removefile(self):
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        
        if len(selected_indexes) == 0:
            return
        
        selectedRows = [x.row() for x in selected_indexes]
        selectedRows.sort(reverse=True)
        
        print(selectedRows)
        
        for row in selectedRows:
            del self.dataClasses[row]
            self.lstSpectra.takeItem(row)
            self.lblStatus.setText("Files removed.")  
        
        
    def savemju(self):
        
        if self.current < 0:
            return
            
        if self.chkSaveinOneFile.isChecked(): #save in one file
        
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                if self.dataClasses[selected_indexes[i].row()].mju != []:
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].energy):
                    l = len(self.dataClasses[y].energy)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if self.dataClasses[y].mju != []:
                    if l>len(self.dataClasses[y].energy):
                        add_array = zeros(l-len(self.dataClasses[y].energy))
                    else:
                        add_array = []
                    save_array.append(concatenate((self.dataClasses[y].energy, add_array)))
                    save_array.append(concatenate((self.dataClasses[y].mju, add_array)))
    

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files

            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return             

            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".mju"
                save_array = []
                save_array.append(self.dataClasses[y].energy)
                save_array.append(self.dataClasses[y].mju)
                savetxt(filename, transpose(save_array))
            
        
    def savexanes(self):
        
        if self.current < 0:
            return
            
        if self.chkSaveinOneFile.isChecked(): #save in one file
        
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                if self.mju[selected_indexes[i].row()] != []:
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].energy):
                    l = len(self.dataClasses[y].energy)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if self.dataClasses[y].mju != []:
                    if l>len(self.dataClasses[y].energy):
                        add_array = zeros(l-len(self.dataClasses[y].energy))
                    else:
                        add_array = []
                    idx = argmin(abs(self.dataClasses[y].energy - self.dataClasses[y].normalizationEnergy))
                    save_array.append(concatenate((self.dataClasses[y].energy, add_array)))
                    save_array.append(concatenate(
                            (self.dataClasses[y].mjuMinusVictoreen / self.dataClasses[y].mju0[idx], add_array)))

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files

            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return 
            
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".xanes"
                idx = argmin(abs(self.dataClasses[y].energy - self.dataClasses[y].normalizationEnergy))
                save_array = []
                save_array.append(self.dataClasses[y].energy)
                save_array.append(self.dataClasses[y].mjuMinusVictoreen / self.dataClasses[y].mju0[idx])
                savetxt(filename, transpose(save_array))

    def saveexafs(self):
        
        if self.current < 0:
            return
            
        if self.chkSaveinOneFile.isChecked(): #save in one file
        
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                column_captions = column_captions + str(selected_items[i].text()) + " "
                column_captions = column_captions + str(selected_items[i].text()) + " "
                
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].k):
                    l = len(self.dataClasses[y].k)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l>len(self.dataClasses[y].k):
                    add_array = zeros(l-len(self.dataClasses[y].k))
                else:
                    add_array = []
                save_array.append(concatenate((self.dataClasses[y].k, add_array)))
                save_array.append(concatenate((self.dataClasses[y].exafsZLC, add_array)))

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files

            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return 
            
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".exafs"
                save_array = []
                save_array.append(self.dataClasses[y].k)
                save_array.append(self.dataClasses[y].exafsZLC)
                savetxt(filename, transpose(save_array))
        

    def saveft(self):
        
        if self.current < 0:
            return
        
        if self.chkSaveinOneFile.isChecked(): #save in one file
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                column_captions = column_captions + str(selected_items[i].text()) + " "
                column_captions = column_captions + str(selected_items[i].text()) + " "
                column_captions = column_captions + str(selected_items[i].text()) + " "
                
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].rsm):
                    l = len(self.dataClasses[y].rsm)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l>len(self.dataClasses[y].rsm):
                    add_array = zeros(l-len(self.dataClasses[y].rsm))
                else:
                    add_array = []
                save_array.append(concatenate((self.dataClasses[y].rsm, add_array)))
                save_array.append(concatenate((self.dataClasses[y].efrsm, add_array)))
                save_array.append(concatenate((self.dataClasses[y].efism, add_array)))

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
        
        else: #save in separate files

            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return 
            
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".ft"
                save_array = []
                save_array.append(self.dataClasses[y].rsm)
                save_array.append(self.dataClasses[y].efrsm)
                save_array.append(self.dataClasses[y].efism)
                savetxt(filename, transpose(save_array))

    def savebft(self):
        
        if self.current < 0:
            return
        
        if self.chkSaveinOneFile.isChecked(): #save in one file
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                column_captions = column_captions + str(selected_items[i].text()) + " "
                column_captions = column_captions + str(selected_items[i].text()) + " "
                
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].bftk):
                    l = len(self.dataClasses[y].bftk)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l>len(self.dataClasses[y].bftk):
                    add_array = zeros(l-len(self.dataClasses[y].bftk))
                else:
                    add_array = []
                save_array.append(concatenate((self.dataClasses[y].bftk, add_array)))
                save_array.append(concatenate((self.dataClasses[y].bftexafs, add_array)))

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files

            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return 
            
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".bft"
                save_array = []
                save_array.append(self.dataClasses[y].bftk)
                save_array.append(self.dataClasses[y].bftexafs)
                savetxt(filename, transpose(save_array))
                
    def saveXes(self):
        if self.current < 0:
            return
            
        if self.chkSaveinOneFile.isChecked(): #save in one file
        
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                if self.mju[selected_indexes[i].row()] != []:
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    column_captions = column_captions + str(selected_items[i].text()) + " "
                    
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].energy):
                    l = len(self.dataClasses[y].energy)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if self.dataClasses[y].xes != []:
                    if l>len(self.dataClasses[y].energy):
                        add_array = zeros(l-len(self.dataClasses[y].energy))
                    else:
                        add_array = []
                    save_array.append(concatenate((self.dataClasses[y].energy, add_array)))
                    save_array.append(concatenate((self.dataClasses[y].xes, add_array)))
    

            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files


            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return             

            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".xes"
                save_array = []
                save_array.append(self.dataClasses[y].energy)
                save_array.append(self.dataClasses[y].xes)
                savetxt(filename, transpose(save_array))
    
    def saveXesAnorm(self):
        if self.current < 0:
            return
            
        if self.chkSaveinOneFile.isChecked(): #save in one file
        
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            column_captions = ""
            for i in range(len(list(selected_items))):
                column_captions = column_captions + str(selected_items[i].text()) + " "
                column_captions = column_captions + str(selected_items[i].text()) + " "
                
            l = 0
            #find the longest array
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l<len(self.dataClasses[y].energy):
                    l = len(self.dataClasses[y].energy)
            
            save_array = []
            for i in range(len(list(selected_indexes))):
                y= selected_indexes[i].row()
                if l>len(self.dataClasses[y].energy):
                    add_array = zeros(l-len(self.dataClasses[y].energy))
                else:
                    add_array = []
                save_array.append(concatenate((self.dataClasses[y].energy, add_array)))
                save_array.append(concatenate((self.dataClasses[y].xesAreaNorm, add_array)))


            fn = self.savefiledialog_qtgui()
            if fn == "":
                return 
            
            savetxt(fn, transpose(save_array), header=column_captions)
            
        else: #save in separate files-
            directory = self.directoryfiledialog_qtgui()
            if directory == "":
                return 
            
            selected_items = self.lstSpectra.selectedItems()
            selected_indexes = self.lstSpectra.selectedIndexes()
    
            for i in range(len(list(selected_items))):
                y= selected_indexes[i].row()
                filename = directory +"/" + str(selected_items[i].text()) + ".xesan"
                save_array = []
                save_array.append(self.dataClasses[y].energy)
                save_array.append(self.dataClasses[y].xesAreaNorm)
                savetxt(filename, transpose(save_array))
                
    def ampPhaSave(self):
        if self.current < 0:
            return
        cnr = self.current
        d = self.directoryfiledialog_qtgui()
        if d == "":
            return 
            
        savetxt(d + "/" + self.lstSpectra.item(cnr).text() + ".amp", 
                transpose([self.dataClasses[cnr].bftk, self.dataClasses[cnr].bftAmp / self.dataClasses[cnr].bftk**self.dataClasses[cnr].kPower ])) 
        
        savetxt(d + "/" + self.lstSpectra.item(cnr).text() + ".pha", 
                transpose([self.dataClasses[cnr].bftk, self.dataClasses[cnr].bftPha ]))

    def lstSpectraItemClicked(self):
        
#        startTime = timer()

        self.current =  self.lstSpectra.currentRow()
        print("Selected element", self.current)
        cnr = self.current
        
        self.lblCurrentFile.setText(self.dataClasses[cnr].name)
        
        if isinstance(self.dataClasses[cnr], xaesa_xes_class):
            #setup interface
            self.frameXasFig.hide()
            self.frameXesFig.show()
            
            self.frameXesParams.show()
            self.frameLParams.hide()
            self.frameLFTParams.hide()
            self.frameLBFTParams.hide()
            
            self.btnCompareXas.setMenu(self.mnuCompareXES)
            self.btnSaveXas.setMenu(self.mnuSaveXES)
            self.btnAverage.setMenu(self.mnuAverageXES)
            self.btnRebin.setMenu(self.mnuRebinXES)
            
            #set params to edit boxes
            self.edtXesE0.setText("{:.2f}".format(self.dataClasses[cnr].E0))
            self.edtXesE1.setText("{:.2f}".format(self.dataClasses[cnr].E1))
            self.edtXesE2.setText("{:.2f}".format(self.dataClasses[cnr].E2))
            self.edtXesE3.setText("{:.2f}".format(self.dataClasses[cnr].E3))
            
            self.edtXesEANormMin.setText("{:.2f}".format(self.dataClasses[cnr].eAreaNormMin))
            self.edtXesEANormMax.setText("{:.2f}".format(self.dataClasses[cnr].eAreaNormMax))
            
            #update graph data
            self.lineXesOriginal.set_xdata(self.dataClasses[cnr].energy)
            self.lineXesOriginal.set_ydata(self.dataClasses[cnr].xes)
            self.lineE0xes.set_xdata(self.dataClasses[cnr].E0)
            self.lineE1xes.set_xdata(self.dataClasses[cnr].E1)
            self.lineE2xes.set_xdata(self.dataClasses[cnr].E2)
            self.lineE3xes.set_xdata(self.dataClasses[cnr].E3)
            self.lineEminAnorm.set_xdata(self.dataClasses[cnr].eAreaNormMin)
            self.lineEmaxAnorm.set_xdata(self.dataClasses[cnr].eAreaNormMax)
            
            self.lineXesBkgrCorr.set_xdata(self.dataClasses[cnr].energy)
            self.lineXesBkgrCorr.set_ydata(self.dataClasses[cnr].xesBkgrCorrected)
            
            self.lineXesNorm.set_xdata(self.dataClasses[cnr].energy)
            self.lineXesNorm.set_ydata(self.dataClasses[cnr].xesAreaNorm)
            
            self.ax_xes.relim()
            self.ax_xes_bc.relim()
            self.ax_xes_norm.relim()
            
            self.ax_xes.autoscale()
            self.ax_xes_bc.autoscale()
            self.ax_xes_norm.autoscale()
            
            self.figXes.tight_layout()
            self.canvXes.draw()
        
        if isinstance(self.dataClasses[cnr], xaesa_exafs_class):
            #setup interface
            self.frameXasFig.show()
            self.frameXesFig.hide()
            
            self.frameXesParams.hide()
            self.frameLParams.show()
            self.frameLFTParams.show()
            self.frameLBFTParams.show()
            
            self.btnCompareXas.setMenu(self.mnuCompareXAS)
            self.btnSaveXas.setMenu(self.mnuSaveXAS)
            self.btnAverage.setMenu(self.mnuAverageXAS)
            self.btnRebin.setMenu(self.mnuRebinXAS)
            
            self.edtE0.setText("{:.2f}".format(self.dataClasses[cnr].E0))
            self.edtE1.setText("{:.2f}".format(self.dataClasses[cnr].E1))
            self.edtE2.setText("{:.2f}".format(self.dataClasses[cnr].E2))
            self.edtE3.setText("{:.2f}".format(self.dataClasses[cnr].E3))
            self.edtkpow.setText(str(self.dataClasses[cnr].kPower))
            self.edtsm.setText(str(self.dataClasses[cnr].zeroLineCorr))
            self.edtmju0poldegree.setText("{:.0f}".format(self.dataClasses[cnr].mju0PolinomialDegree))
    
            self.edtkmin.setText(str(self.dataClasses[cnr].kMin))
            self.edtkmax.setText("{:.2f}".format(self.dataClasses[cnr].kMax))
            self.edtdk.setText(str(self.dataClasses[cnr].dk))
            self.edtrmin.setText(str(self.dataClasses[cnr].rMin))
            self.edtrmax.setText(str(self.dataClasses[cnr].rMax))
            self.edtdr.setText(str(self.dataClasses[cnr].dr))
    
            self.edtrminbft.setText(str(self.dataClasses[cnr].rMinBft))
            self.edtrmaxbft.setText(str(self.dataClasses[cnr].rMaxBft))
            self.edtbftwindowparam.setText(str(self.dataClasses[cnr].bftWindowParam))
            
            self.edtExafsNormEnergy.setText(str(self.dataClasses[cnr].normalizationEnergy))
            self.chkExafsNormalization.setChecked(self.dataClasses[cnr].normalizationMode)
            
            self.abs_scatter.set_xdata(self.dataClasses[cnr].energy)
            self.abs_scatter.set_ydata(self.dataClasses[cnr].mju)
            self.mjub_line.set_xdata(self.dataClasses[cnr].energy)
            self.mjub_line.set_ydata(self.dataClasses[cnr].victoreen)
            self.mju0_mjub_line.set_xdata(self.dataClasses[cnr].energy)
            self.mju0_mjub_line.set_ydata(self.dataClasses[cnr].mju0+self.dataClasses[cnr].victoreen)
            self.deriv_line.set_xdata(self.dataClasses[cnr].energy)
            self.deriv_line.set_ydata(self.dataClasses[cnr].mjuDerivative)
            self.lineE0.set_xdata(self.dataClasses[cnr].E0)
            self.lineE1.set_xdata(self.dataClasses[cnr].E1)
            self.lineE2.set_xdata(self.dataClasses[cnr].E2)
            self.lineE3.set_xdata(self.dataClasses[cnr].E3)
            self.ax_abs.relim()
            self.ax_abs.autoscale()
            self.ax_abs2.relim()
            self.ax_abs2.autoscale()            
    
            self.exafs_line.set_xdata(self.dataClasses[cnr].k)  
            self.exafs_line.set_ydata(self.dataClasses[cnr].exafs)
            self.exafssm_line.set_xdata(self.dataClasses[cnr].k)
            self.exafssm_line.set_ydata(self.dataClasses[cnr].exafsZeroLine)
            self.ax_exafs.relim()
            self.ax_exafs.autoscale()
    
            self.efr_line.set_xdata(self.dataClasses[cnr].r)
            self.efr_line.set_ydata(self.dataClasses[cnr].efr)
            self.efi_line.set_xdata(self.dataClasses[cnr].r)
            self.efi_line.set_ydata(self.dataClasses[cnr].efi)
            
            self.efrsm_line.set_xdata(self.dataClasses[cnr].rZLC)
            self.efrsm_line.set_ydata(self.dataClasses[cnr].efrZLC)
            self.efism_line.set_xdata(self.dataClasses[cnr].rZLC)
            self.efism_line.set_ydata(self.dataClasses[cnr].efiZLC)
            
            self.bftwinr_line.set_xdata(self.dataClasses[cnr].rZLC)
            self.bftwinr_line.set_ydata(self.dataClasses[cnr].bftefrWindow)
            self.bftwini_line.set_xdata(self.dataClasses[cnr].rZLC)
            self.bftwini_line.set_ydata(self.dataClasses[cnr].bftefiWindow)
            
            self.ax_ft.relim()
            self.ax_ft.autoscale()
    
            self.exafsd_line.set_xdata(self.dataClasses[cnr].k)
            self.exafsd_line.set_ydata(self.dataClasses[cnr].exafsZLC)
            
            self.bftexafs_line.set_xdata(self.dataClasses[cnr].bftk)
            self.bftexafs_line.set_ydata(self.dataClasses[cnr].bftEXAFS)
            
            self.ax_exafssm.relim()
            self.ax_exafssm.autoscale()
            self.fig.tight_layout()
            self.canv.draw()
            
    #        endTime = timer()
    #        self.lblStatus.setText("Replot time {:.4f} seconds".format(endTime - startTime) )
            self.lblnpoints.setText("Number of experimental points : " + str(len(self.dataClasses[cnr].mju)))
        
        
    def listItemRightClicked(self, QPos): 
        
        parentPosition = self.lstSpectra.mapToGlobal(QtCore.QPoint(0, 0))        
        self.lstSpectraMenu.move(parentPosition + QPos)
        self.lstSpectraMenu.show() 
            
    def lstSpectraItemDoubleClicked(self):
        item = self.lstSpectra.currentItem()
        item.setFlags(item.flags() | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsEnabled)
        self.lstSpectra.editItem(self.lstSpectra.currentItem())

    def savehdf5(self):

        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setAcceptMode(1) # save dialog
        dlg.setNameFilters(["HDF5 files (*.hdf5)",  "All files (*.*)"])
#        dlg.setDirectory(self.currenthdf5)
#        dlg.setDirectory(self.currentdir)
        if dlg.exec_():
            fnlist = dlg.selectedFiles()
        else:
            return
        fn = fnlist[0]    
        head, tail = path.split(str(fn))
#        self.currenthdf5 = head
        self.currentdir = head

        fhdf5 = h5py.File(fn, "w")

#        anfn = []
#        for i in range(0, self.lstSpectra.count()):
#            anfn.append(str(self.lstSpectra.item(i).text()))
#        asciiList = [n.encode("ascii", "ignore") for n in anfn]

                     
        for i in range(self.lstSpectra.count()):
            
            if isinstance(self.dataClasses[i], xaesa_exafs_class): #save XAS data                
            
                grp = fhdf5.create_group("{0:06d}".format(i) + "&" + self.lstSpectra.item(i).text())
                grp.create_dataset("classType", data = 0)
                
                grp.create_dataset("raw_data_type", data = self.dataClasses[i].raw_data_type)
                
                grp.create_dataset("name", data = self.dataClasses[i].name)
                grp.create_dataset("energy", data = self.dataClasses[i].energy)
                grp.create_dataset("energyRebined", data = self.dataClasses[i].energyRebined)
                grp.create_dataset("energyOriginal", data = self.dataClasses[i].energyOriginal)
                
                grp.create_dataset("i0", data = self.dataClasses[i].i0)
                grp.create_dataset("i1", data = self.dataClasses[i].i1)
                grp.create_dataset("i2", data = self.dataClasses[i].i2)
                grp.create_dataset("ifluo", data = self.dataClasses[i].ifluo)
                
                #Mju datasets
                grp.create_dataset("mju", data = self.dataClasses[i].mju)
                grp.create_dataset("mjuRebined", data = self.dataClasses[i].mjuRebined)
                grp.create_dataset("mjuOriginal", data = self.dataClasses[i].mjuOriginal)
                grp.create_dataset("mjuDerivative", data = self.dataClasses[i].mjuDerivative)
                grp.create_dataset("victoreen", data = self.dataClasses[i].victoreen)
                grp.create_dataset("mjuMinusVictoreen", data = self.dataClasses[i].mjuMinusVictoreen)
                grp.create_dataset("mju0", data = self.dataClasses[i].mju0)       
                
                #EXAFS datasets
                grp.create_dataset("k", data = self.dataClasses[i].k )
                grp.create_dataset("exafs", data = self.dataClasses[i].exafs )
                grp.create_dataset("exafsZeroLine", data = self.dataClasses[i].exafsZeroLine )
                grp.create_dataset("exafsZLC", data = self.dataClasses[i].exafsZLC )
                
                # FT datasets
                grp.create_dataset("window", data = self.dataClasses[i].window )
                grp.create_dataset("exafsTimesWindow", data = self.dataClasses[i].exafsTimesWindow )
                grp.create_dataset("exafsZLCTimesWindow", data = self.dataClasses[i].exafsZLCTimesWindow )
                grp.create_dataset("r", data = self.dataClasses[i].r )
                grp.create_dataset("fr", data = self.dataClasses[i].fr )
                grp.create_dataset("fi", data = self.dataClasses[i].fi )       
                grp.create_dataset("efr", data = self.dataClasses[i].efr )
                grp.create_dataset("efi", data = self.dataClasses[i].efi )
                
                grp.create_dataset("rZLC", data = self.dataClasses[i].rZLC )
                grp.create_dataset("frZLC", data = self.dataClasses[i].frZLC )
                grp.create_dataset("fiZLC", data = self.dataClasses[i].fiZLC )       
                grp.create_dataset("efrZLC", data = self.dataClasses[i].efrZLC )
                grp.create_dataset("efiZLC", data = self.dataClasses[i].efiZLC )
                
                # BFT datasets
                grp.create_dataset("bftWindow", data = self.dataClasses[i].bftWindow )
                grp.create_dataset("bftk", data = self.dataClasses[i].bftk )
                grp.create_dataset("bftefr", data = self.dataClasses[i].bftefr )
                grp.create_dataset("bftefi", data = self.dataClasses[i].bftefi )
                grp.create_dataset("bftAmp", data = self.dataClasses[i].bftAmp )
                grp.create_dataset("bftPha", data = self.dataClasses[i].bftPha )
                grp.create_dataset("bftEXAFS", data = self.dataClasses[i].bftEXAFS )
                grp.create_dataset("bftefrWindow", data = self.dataClasses[i].bftefrWindow )
                grp.create_dataset("bftefiWindow", data = self.dataClasses[i].bftefiWindow )    
                
                #Rebin parameters
                grp.create_dataset("dE1", data = self.dataClasses[i].dE1 )
                grp.create_dataset("dE2", data = self.dataClasses[i].dE2 )
                grp.create_dataset("dE3", data = self.dataClasses[i].dE3 )
                
                
                #extraction params
                grp.create_dataset("E0", data = self.dataClasses[i].E0 )
                grp.create_dataset("E1", data = self.dataClasses[i].E1 )
                grp.create_dataset("E2", data = self.dataClasses[i].E2 )
                grp.create_dataset("E3", data = self.dataClasses[i].E3 )
                
                grp.create_dataset("kPower", data = self.dataClasses[i].kPower)
                
                grp.create_dataset("zeroLineCorr", data = self.dataClasses[i].zeroLineCorr )
                
                grp.create_dataset("mju0PolinomialDegree", data = self.dataClasses[i].mju0PolinomialDegree )
                
                grp.create_dataset("normalizationMode", data = self.dataClasses[i].normalizationMode ) #0 for mju0 normalization, 1 for value normalization at given energy
                grp.create_dataset("normalizationEnergy", data = self.dataClasses[i].normalizationEnergy )
                
                #FT params
                grp.create_dataset("kMin", data = self.dataClasses[i].kMin )
                grp.create_dataset("kMax", data = self.dataClasses[i].kMax )
                grp.create_dataset("dk", data = self.dataClasses[i].dk )
                grp.create_dataset("rMin", data = self.dataClasses[i].rMin )
                grp.create_dataset("rMax", data = self.dataClasses[i].rMax )
                grp.create_dataset("dr", data = self.dataClasses[i].dr )
                
                #BFT params
                grp.create_dataset("rMinBft", data = self.dataClasses[i].rMinBft )
                grp.create_dataset("rMaxBft", data = self.dataClasses[i].rMaxBft )
                grp.create_dataset("bftWindowParam", data = self.dataClasses[i].bftWindowParam )       
                
                # Deglitching params and result
                
                #fit params and results
                
                grp.create_dataset("fitKMin", data = self.dataClasses[i].fitKMin )
                grp.create_dataset("fitKMax", data = self.dataClasses[i].fitKMax )
                grp.create_dataset("fitdk", data = self.dataClasses[i].fitdk )
                
                grp.create_dataset("fitNShels", data = self.dataClasses[i].fitNShels )
                
                grp.create_dataset("fitParams", data = self.dataClasses[i].fitParams )
                grp.create_dataset("fitAmps", data = self.dataClasses[i].fitAmps )
                grp.create_dataset("fitPhas", data = self.dataClasses[i].fitPhas )
                
                grp.create_dataset("fitK", data = self.dataClasses[i].fitK )
                grp.create_dataset("fitExafs", data = self.dataClasses[i].fitExafs )
                
                #rdf params and results
                
                grp.create_dataset("isRdfed", data = self.dataClasses[i].isRdfed )
                
                grp.create_dataset("rdfKMin", data = self.dataClasses[i].rdfKMin )
                grp.create_dataset("rdfKMax", data = self.dataClasses[i].rdfKMax )
                grp.create_dataset("rdfdk", data = self.dataClasses[i].rdfdk )
                
                grp.create_dataset("rdfRMin", data = self.dataClasses[i].rdfRMin )
                grp.create_dataset("rdfRMax", data = self.dataClasses[i].rdfRMax )
                grp.create_dataset("rdfdr", data = self.dataClasses[i].rdfdr )
                
                grp.create_dataset("rdfMaxIterations", data = self.dataClasses[i].rdfMaxIterations )
                
                grp.create_dataset("rdfAmpK", data = self.dataClasses[i].rdfAmpK )
                grp.create_dataset("rdfAmp", data = self.dataClasses[i].rdfAmp )
                grp.create_dataset("rdfPhaK", data = self.dataClasses[i].rdfPhaK )
                grp.create_dataset("rdfPha", data = self.dataClasses[i].rdfPha )
                
                grp.create_dataset("rdfAmpFile", data = self.dataClasses[i].rdfAmpFile )
                grp.create_dataset("rdfPhaFile", data = self.dataClasses[i].rdfPhaFile )
                
                grp.create_dataset("rdfK", data = self.dataClasses[i].rdfK )
                grp.create_dataset("rdfExafs", data = self.dataClasses[i].rdfExafs )
                
                grp.create_dataset("rdfR", data = self.dataClasses[i].rdfR )
                grp.create_dataset("rdf", data = self.dataClasses[i].rdf )
                
            
            if isinstance(self.dataClasses[i], xaesa_xes_class): #save XES data                
            
                grp = fhdf5.create_group("{0:06d}".format(i) + "&" + self.lstSpectra.item(i).text())
                grp.create_dataset("classType", data = 10)
                grp.create_dataset("name", data = self.dataClasses[i].name)
                grp.create_dataset("energy", data = self.dataClasses[i].energy)
                grp.create_dataset("energyRebined", data = self.dataClasses[i].energyRebined)
                grp.create_dataset("energyOriginal", data = self.dataClasses[i].energyOriginal)
                
                grp.create_dataset("xes", data = self.dataClasses[i].xes)
                grp.create_dataset("xesBackground", data = self.dataClasses[i].xesBackground)
                grp.create_dataset("xesBkgrCorrected", data = self.dataClasses[i].xesBkgrCorrected)
                grp.create_dataset("xesAreaNorm", data = self.dataClasses[i].xesAreaNorm)
                grp.create_dataset("xesMaxNorm", data = self.dataClasses[i].xesMaxNorm)
                grp.create_dataset("xesRebinned", data = self.dataClasses[i].xesRebinned)
                grp.create_dataset("xesOriginal", data = self.dataClasses[i].xesOriginal)
                
                grp.create_dataset("E0", data = self.dataClasses[i].E0)
                grp.create_dataset("E1", data = self.dataClasses[i].E1)
                grp.create_dataset("E2", data = self.dataClasses[i].E2)
                grp.create_dataset("E3", data = self.dataClasses[i].E3)
                
                grp.create_dataset("eAreaNormMin", data = self.dataClasses[i].eAreaNormMin)
                grp.create_dataset("eAreaNormMax", data = self.dataClasses[i].eAreaNormMax)                

        fhdf5.close()

        return

    def openhdf5(self):
            
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFile)
        dlg.setAcceptMode(0) # open dialog
        dlg.setNameFilters(["HDF5 files (*.hdf5)",  "All files (*.*)"])
#        dlg.setDirectory(self.currenthdf5)
        dlg.setDirectory(self.currentdir)
        if dlg.exec_():
            fnlist = dlg.selectedFiles()
        else:
            return
        fn = fnlist[0]
        head, tail = path.split(str(fn))
            
#        self.currenthdf5 = head
        self.currentdir = head

        f = h5py.File(fn, "r")
        
        self.lstSpectra.clear()
        self.dataClasses = []
        
        
        mainKeyList = list(f.keys())
        
        for i in range(0, len(mainKeyList)):
            #remove numbering if existing
            if '&' in mainKeyList[i]:
                txt = mainKeyList[i].split('&', 1)[-1]
                self.lstSpectra.addItem(txt)
            else:
                self.lstSpectra.addItem(mainKeyList[i])

            grp = f.get(mainKeyList[i])
            
            classType = grp.get("classType").value
            print("classType", classType)
        ##################### OPEN XAS    
            if classType == 0: #XAS
                
                
                
                self.dataClasses.append(xaesa_exafs_class(grp.get("raw_data_type").value))
                self.dataClasses[-1].raw_data_type = grp.get("raw_data_type").value
                self.dataClasses[-1].name = grp.get("name").value
                self.dataClasses[-1].energy = grp.get("energy").value
                self.dataClasses[-1].energyRebined = grp.get("energyRebined").value
                self.dataClasses[-1].energyOriginal = grp.get("energyOriginal").value
                
                self.dataClasses[-1].i0 = grp.get("i0").value
                self.dataClasses[-1].i1 = grp.get("i1").value
                self.dataClasses[-1].i2 = grp.get("i2").value
                self.dataClasses[-1].ifluo = grp.get("ifluo").value
                
                
                #Mju datasets
                self.dataClasses[-1].mju = grp.get("mju").value
                self.dataClasses[-1].mjuRebined = grp.get("mjuRebined").value
                self.dataClasses[-1].mjuOriginal = grp.get("mjuOriginal").value
                self.dataClasses[-1].mjuDerivative= grp.get("mjuDerivative").value
                self.dataClasses[-1].victoreen = grp.get("victoreen").value
                self.dataClasses[-1].mjuMinusVictoreen = grp.get("mjuMinusVictoreen").value
                self.dataClasses[-1].mju0 = grp.get("mju0").value      
                
                #EXAFS datasets
                self.dataClasses[-1].k = grp.get("k").value
                self.dataClasses[-1].exafs = grp.get("exafs").value
                self.dataClasses[-1].exafsZeroLine = grp.get("exafsZeroLine").value
                self.dataClasses[-1].exafsZLC = grp.get("exafsZLC").value
                
                # FT datasets
                self.dataClasses[-1].window = grp.get("window").value
                self.dataClasses[-1].exafsTimesWindow = grp.get("exafsTimesWindow").value
                self.dataClasses[-1].exafsZLCTimesWindow = grp.get("exafsZLCTimesWindow").value
                self.dataClasses[-1].r = grp.get("r").value
                self.dataClasses[-1].fr = grp.get("fr").value
                self.dataClasses[-1].fi = grp.get("fi").value        
                self.dataClasses[-1].efr = grp.get("efr").value
                self.dataClasses[-1].efi = grp.get("efi").value
                
                self.dataClasses[-1].rZLC = grp.get("rZLC").value
                self.dataClasses[-1].frZLC = grp.get("frZLC").value
                self.dataClasses[-1].fiZLC = grp.get("fiZLC").value       
                self.dataClasses[-1].efrZLC = grp.get("efrZLC").value
                self.dataClasses[-1].efiZLC = grp.get("efiZLC").value
                
                # BFT datasets
                self.dataClasses[-1].bftWindow = grp.get("bftWindow").value
                self.dataClasses[-1].bftk = grp.get("bftk").value
                self.dataClasses[-1].bftefr = grp.get("bftefr").value
                self.dataClasses[-1].bftefi = grp.get("bftefi").value
                self.dataClasses[-1].bftAmp = grp.get("bftAmp").value
                self.dataClasses[-1].bftPha = grp.get("bftPha").value
                self.dataClasses[-1].bftEXAFS = grp.get("bftEXAFS").value
                self.dataClasses[-1].bftefrWindow = grp.get("bftefrWindow").value
                self.dataClasses[-1].bftefiWindow = grp.get("bftefiWindow").value    
                
                #Rebin parameters
                self.dataClasses[-1].dE1 = grp.get("dE1").value
                self.dataClasses[-1].dE2 = grp.get("dE2").value
                self.dataClasses[-1].dE3 = grp.get("dE3").value
                
                
                #extraction params
                self.dataClasses[-1].E0 = grp.get("E0").value
                self.dataClasses[-1].E1 = grp.get("E1").value
                self.dataClasses[-1].E2 = grp.get("E2").value
                self.dataClasses[-1].E3 = grp.get("E3").value
                
                self.dataClasses[-1].kPower = grp.get("kPower").value
                
                self.dataClasses[-1].zeroLineCorr = grp.get("zeroLineCorr").value
                
                self.dataClasses[-1].mju0PolinomialDegree = grp.get("mju0PolinomialDegree").value
                
                self.dataClasses[-1].normalizationMode = grp.get("normalizationMode").value #0 for mju0 normalization, 1 for value normalization at given energy
                self.dataClasses[-1].normalizationEnergy = grp.get("normalizationEnergy").value
                
                #FT params
                self.dataClasses[-1].kMin = grp.get("kMin").value
                self.dataClasses[-1].kMax = grp.get("kMax").value
                self.dataClasses[-1].dk = grp.get("dk").value
                self.dataClasses[-1].rMin = grp.get("rMin").value
                self.dataClasses[-1].rMax = grp.get("rMax").value
                self.dataClasses[-1].dr = grp.get("dr").value
                
                #BFT params
                self.dataClasses[-1].rMinBft = grp.get("rMinBft").value
                self.dataClasses[-1].rMaxBft = grp.get("rMaxBft").value
                self.dataClasses[-1].bftWindowParam = grp.get("bftWindowParam").value       
                
                # Deglitching params and result
                
                #fit params and results
                
                self.dataClasses[-1].fitKMin = grp.get("fitKMin").value
                self.dataClasses[-1].fitKMax = grp.get("fitKMax").value
                self.dataClasses[-1].fitdk = grp.get("fitdk").value
                
                self.dataClasses[-1].fitNShels = grp.get("fitNShels").value
                
                self.dataClasses[-1].fitParams = grp.get("fitParams").value
                self.dataClasses[-1].fitAmps = grp.get("fitAmps").value
                self.dataClasses[-1].fitPhas = grp.get("fitPhas").value
                
                self.dataClasses[-1].fitK = grp.get("fitK").value
                self.dataClasses[-1].fitExafs = grp.get("fitExafs").value
                
                #rdf params and results
                
                self.dataClasses[-1].isRdfed = grp.get("isRdfed").value
                
                self.dataClasses[-1].rdfKMin = grp.get("rdfKMin").value
                self.dataClasses[-1].rdfKMax = grp.get("rdfKMax").value
                self.dataClasses[-1].rdfdk = grp.get("rdfdk").value
                
                self.dataClasses[-1].rdfRMin = grp.get("rdfRMin").value
                self.dataClasses[-1].rdfRMax =grp.get("rdfRMax").value
                self.dataClasses[-1].rdfdr = grp.get("rdfdr").value
                
                self.dataClasses[-1].rdfMaxIterations = grp.get("rdfMaxIterations").value
                
                self.dataClasses[-1].rdfAmpK = grp.get("rdfAmpK").value
                self.dataClasses[-1].rdfAmp = grp.get("rdfAmp").value
                self.dataClasses[-1].rdfPhaK = grp.get("rdfPhaK").value
                self.dataClasses[-1].rdfPha = grp.get("rdfPha").value
                
                self.dataClasses[-1].rdfAmpFile = grp.get("rdfAmpFile").value
                self.dataClasses[-1].rdfPhaFile = grp.get("rdfPhaFile").value
                
                self.dataClasses[-1].rdfK = grp.get("rdfK").value
                self.dataClasses[-1].rdfExafs = grp.get("rdfExafs").value
                
                self.dataClasses[-1].rdfR = grp.get("rdfR").value
                self.dataClasses[-1].rdf = grp.get("rdf").value
                
                sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                
                if grp.get("raw_data_type").value == 3:
                    sel_item.setForeground(QtCore.Qt.blue)
                else:
                    sel_item.setForeground(QtCore.Qt.darkGreen)
                
        ##################### OPEN XES
            
            if classType == 10: #XES
                self.dataClasses.append(xaesa_xes_class())
                self.dataClasses[-1].name = grp.get("name").value
                self.dataClasses[-1].energy = grp.get("energy").value
                self.dataClasses[-1].energyRebined = grp.get("energyRebined").value
                self.dataClasses[-1].energyOriginal = grp.get("energyOriginal").value
                
                self.dataClasses[-1].xes = grp.get("xes").value
                self.dataClasses[-1].xesBackground = grp.get("xesBackground").value
                self.dataClasses[-1].xesBkgrCorrected = grp.get("xesBkgrCorrected").value
                self.dataClasses[-1].xesAreaNorm = grp.get("xesAreaNorm").value
                self.dataClasses[-1].xesMaxNorm = grp.get("xesMaxNorm").value
                self.dataClasses[-1].xesRebinned = grp.get("xesRebinned").value
                self.dataClasses[-1].xesOriginal = grp.get("xesOriginal").value
                
                self.dataClasses[-1].E0 = grp.get("E0").value
                self.dataClasses[-1].E1 = grp.get("E1").value
                self.dataClasses[-1].E2 = grp.get("E2").value
                self.dataClasses[-1].E3 = grp.get("E3").value
                self.dataClasses[-1].eAreaNormMin = grp.get("eAreaNormMin").value
                self.dataClasses[-1].eAreaNormMax = grp.get("eAreaNormMax").value
                
                sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
                sel_item.setForeground(QtCore.Qt.darkMagenta)
  
   
        f.close()

#        self.current = 0
#        self.lstSpectra.setCurrentRow(0)        
#        self.lstSpectraItemClicked()

    def CopyParams(self):
        cnr = self.current
        
        if isinstance(self.dataClasses[cnr], xaesa_xes_class):
            self.copiedparamsXES = [self.dataClasses[cnr].E0,
                                 self.dataClasses[cnr].E1,
                                 self.dataClasses[cnr].E2,
                                 self.dataClasses[cnr].E3,
                                 self.dataClasses[cnr].eAreaNormMin,
                                 self.dataClasses[cnr].eAreaNormMax
                                 ]
            print("XES params copied")
            
        if isinstance(self.dataClasses[cnr], xaesa_exafs_class):    
        
            self.copiedparamsXAS = [self.dataClasses[cnr].E0,
                                 self.dataClasses[cnr].E1,
                                 self.dataClasses[cnr].E2,
                                 self.dataClasses[cnr].E3,
                                 self.dataClasses[cnr].zeroLineCorr,
                                 self.dataClasses[cnr].kPower,
                                 self.dataClasses[cnr].mju0PolinomialDegree,
                                 self.dataClasses[cnr].normalizationMode,
                                 self.dataClasses[cnr].normalizationEnergy,
                                 self.dataClasses[cnr].kMin,
                                 self.dataClasses[cnr].kMax,
                                 self.dataClasses[cnr].dk ,
                                 self.dataClasses[cnr].rMin,
                                 self.dataClasses[cnr].rMax,
                                 self.dataClasses[cnr].dr,                    
                                 self.dataClasses[cnr].rMinBft,
                                 self.dataClasses[cnr].rMaxBft,
                                 self.dataClasses[cnr].bftWindowParam]
        
            print("XAS params copied")
    
    def ApplytoSelected(self):
        if len(self.copiedparamsXAS)==0 and len(self.copiedparamsXES)==0:
            return
        selected_indexes = self.lstSpectra.selectedIndexes()
        
        startTime = timer()

        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            print(y)
            
            if isinstance(self.dataClasses[y], xaesa_xes_class):
                self.dataClasses[y].E0 = self.copiedparamsXES[0]
                self.dataClasses[y].E1 = self.copiedparamsXES[1]
                self.dataClasses[y].E2 = self.copiedparamsXES[2]
                self.dataClasses[y].E3 = self.copiedparamsXES[3]
                
                self.dataClasses[y].eAreaNormMin = self.copiedparamsXES[4]
                self.dataClasses[y].eAreaNormMax = self.copiedparamsXES[5]
                
                self.dataClasses[y].removeBackground()
                self.dataClasses[y].areaNormalize()
                
                print("XES params applied")
            
            if isinstance(self.dataClasses[y], xaesa_exafs_class):
                self.dataClasses[y].E0 = self.copiedparamsXAS[0]
                self.dataClasses[y].E1 = self.copiedparamsXAS[1]
                self.dataClasses[y].E2 = self.copiedparamsXAS[2]
                self.dataClasses[y].E3 = self.copiedparamsXAS[3]
                self.dataClasses[y].zeroLineCorr = self.copiedparamsXAS[4]
                self.dataClasses[y].kPower = self.copiedparamsXAS[5]
                self.dataClasses[y].mju0PolinomialDegree = self.copiedparamsXAS[6]
    
                self.dataClasses[y].normalizationMode = self.copiedparamsXAS[7]
                self.dataClasses[y].normalizationEnergy = self.copiedparamsXAS[8]
        
                self.dataClasses[y].kMin = self.copiedparamsXAS[9]
                self.dataClasses[y].kMax = self.copiedparamsXAS[10]
                self.dataClasses[y].dk = self.copiedparamsXAS[11]
                self.dataClasses[y].rMin = self.copiedparamsXAS[12]
                self.dataClasses[y].rMax = self.copiedparamsXAS[13]
                self.dataClasses[y].dr = self.copiedparamsXAS[14]
        
                self.dataClasses[y].rMinBft = self.copiedparamsXAS[15]
                self.dataClasses[y].rMaxBft = self.copiedparamsXAS[16]
                self.dataClasses[y].bftWindowParam = self.copiedparamsXAS[17]
                
                self.dataClasses[y].redoExtraction()
                
                print("XAS params applied")
            
        endTime = timer()
        self.lblStatus.setText("Params applied to {:d} files in {:.4f} seconds".format(len(list(selected_indexes)), endTime - startTime) )
            
    def removeglitches(self):
        
        if self.current < 0:
            return
        
        cnr = self.current
        
        self.w = DGWindow()
        self.w.exafs = copy(self.dataClasses[cnr].exafsZLC)
        self.w.exafsdg = copy(self.dataClasses[cnr].exafsZLC)
        self.w.k = copy(self.dataClasses[cnr].k)
        self.w.plot()

        if self.w.exec_():
            self.dataClasses[cnr].exafsZLCdeglitch = copy(self.w.exafsdg) 
            self.dataClasses[cnr].exafsZLC = copy(self.w.exafsdg) 
            self.dataClasses[cnr].redoFtBft()
#            self.glitchesRemoved[cnr] = 1
            
#        self.ees_ft()
#        self.lstSpectraItemClicked()
        
    def fit(self):
        
        if self.current < 0:
            return
        
        cnr = self.current
        fw = FitWindow()
        
        if self.dataClasses[cnr].isFitted == 1:
            fw.fit_result = copy(self.dataClasses[cnr].fitExafs)
            fw.savedshellnr = self.dataClasses[cnr].fitNShels
            fw.fit_amps = copy(self.dataClasses[cnr].fitAmps)
            fw.fit_phases = copy(self.dataClasses[cnr].fitPhas)
            fw.fit_params = copy(self.dataClasses[cnr].fitParams)
            fw.costfunction = 0 #self.fit_costfunction[cnr]
            fw.ksettings = [[self.dataClasses[cnr].fitKMin, self.dataClasses[cnr].fitKMax, self.dataClasses[cnr].fitdk]]
            fw.isfitted = 1
            fw.updateUI()
        
        fw.bft = copy(self.dataClasses[cnr].bftEXAFS)
        fw.k = copy(self.dataClasses[cnr].bftk)
        fw.bftft()
        fw.plot()
        
        if fw.exec_(): #copy results and params to arrays
            self.dataClasses[cnr].isFitted = 1
            self.dataClasses[cnr].fitExafs = copy(fw.fit_result)
            self.dataClasses[cnr].fitNShels = fw.shellnr
            self.dataClasses[cnr].fitAmps = copy(fw.fit_amps)
            self.dataClasses[cnr].fitPhas = copy(fw.fit_phases)
            self.dataClasses[cnr].fitParams = copy(fw.fit_params)
#            self.dataClasses[cnr].fit_costfunction[cnr] = fw.costfunction
            self.dataClasses[cnr].fitKMin = fw.ksettings[0][0]
            self.dataClasses[cnr].fitKMax = fw.ksettings[0][1]
            self.dataClasses[cnr].fitdk = fw.ksettings[0][2]
            pass
        
    def rdf(self):
        
        if self.current < 0:
            return
        
        cnr = self.current
        self.rdfw = RdfWindow()
        
        if self.dataClasses[cnr].isRdfed == 1 :
            self.rdfw.rdf = copy(self.dataClasses[cnr].rdf)
            self.rdfw.rdf_r = copy(self.dataClasses[cnr].rdfR)
            self.rdfw.rdf_exafs_k = copy(self.dataClasses[cnr].rdfK)
            self.rdfw.rdf_exafs = copy(self.dataClasses[cnr].rdfExafs)
            self.rdfw.params = [self.dataClasses[cnr].rdfKMin,
                                self.dataClasses[cnr].rdfKMax,
                                self.dataClasses[cnr].rdfdk,
                                self.dataClasses[cnr].rdfRMin,
                                self.dataClasses[cnr].rdfRMax,
                                self.dataClasses[cnr].rdfdr,
                                self.dataClasses[cnr].rdfMaxIterations            
                               ]
            
            self.rdfw.amp_orig[0] = copy(self.dataClasses[cnr].rdfAmpK)
            self.rdfw.amp_orig[1] = copy(self.dataClasses[cnr].rdfAmp) 
            self.rdfw.pha_orig[0] = copy(self.dataClasses[cnr].rdfPhaK)
            self.rdfw.pha_orig[1] = copy(self.dataClasses[cnr].rdfPha)
        
        self.rdfw.bft = copy(self.dataClasses[cnr].bftEXAFS)
        self.rdfw.k = copy(self.dataClasses[cnr].bftk)
        self.rdfw.kpow = self.dataClasses[cnr].kPower
        self.rdfw.bftft()
        self.rdfw.plot()
        
        myStream = MyStream()
        myStream.message.connect(self.rdfw.on_myStream_message)

#        sys.stdout = myStream 
  
        
        if self.rdfw.exec_():
            self.dataClasses[cnr].isRdfed = 1
            self.dataClasses[cnr].rdf = copy(self.rdfw.rdf)
            self.dataClasses[cnr].rdfR = copy(self.rdfw.rdf_r)
            self.dataClasses[cnr].rdfKMin = self.rdfw.params[0]
            self.dataClasses[cnr].rdfKMax = self.rdfw.params[1]
            self.dataClasses[cnr].rdfdk = self.rdfw.params[2]
            self.dataClasses[cnr].rdfRMin = self.rdfw.params[3]
            self.dataClasses[cnr].rdfRax = self.rdfw.params[4]
            self.dataClasses[cnr].rdfdr = self.rdfw.params[5]
            self.dataClasses[cnr].rdfMaxIterations = self.rdfw.params[6]          
            self.dataClasses[cnr].rdfK = copy(self.rdfw.rdf_exafs_k)
            self.dataClasses[cnr].rdfExafs = copy(self.rdfw.rdf_exafs)
            self.dataClasses[cnr].rdfAmpK = copy(self.rdfw.amp_orig[0])
            self.dataClasses[cnr].rdfAmp = copy(self.rdfw.amp_orig[1])
            self.dataClasses[cnr].rdfPhaK = copy(self.rdfw.pha_orig[0])
            self.dataClasses[cnr].rdfPha = copy(self.rdfw.pha_orig[1])
            pass
        
    def compare(self, dataType):
        
        viewerDialog = xaesaViewerWindow()
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return
        
        if  dataType == 'i0' or dataType == 'i1' or dataType == 'i2' or dataType == 'ifluo':
            viewerDialog.viewer.x_caption = 'Energy, eV'
            viewerDialog.viewer.y_caption = 'Intensity'
            
        if  dataType == 'mju' or dataType == 'mjuRebinVSoriginal':
            viewerDialog.viewer.x_caption = 'Energy, eV'
            viewerDialog.viewer.y_caption = 'Absorption, a.u.'
        
        for i in range(len(list(selected_items))):
            
            y= selected_indexes[i].row()
            
            if dataType == 'i0' or dataType == 'i1' or dataType == 'i2':
                if self.dataClasses[y].energyOriginal != []:
                    energy = self.dataClasses[y].energyOriginal
                else:
                    energy = self.dataClasses[y].energy
                viewerDialog.viewer.x_data.append(energy)
                if dataType == 'i0': viewerDialog.viewer.y_data.append(self.dataClasses[y].i0)
                if dataType == 'i1': viewerDialog.viewer.y_data.append(self.dataClasses[y].i1)
                if dataType == 'i2': viewerDialog.viewer.y_data.append(self.dataClasses[y].i2)
                viewerDialog.viewer.labels.append(str(selected_items[i].text()))

                
            if dataType == 'mjuRebinVSoriginal':
                viewerDialog.viewer.x_data.append(self.dataClasses[y].energyOriginal)
                viewerDialog.viewer.y_data.append(self.dataClasses[y].mjuOriginal)
                viewerDialog.viewer.labels.append(str(selected_items[i].text() + 'Original'))
                viewerDialog.viewer.x_data.append(self.dataClasses[y].energy)
                viewerDialog.viewer.y_data.append(self.dataClasses[y].mju)
                viewerDialog.viewer.labels.append(str(selected_items[i].text() + 'Rebin'))

        viewerDialog.viewer.plot()        
        
        viewerDialog.exec_()
        
    def comparemju(self):
        
        self.cw = CompareWindow()
        #self.cw.exafs = np.copy(self.exafsd[cnr])
        #self.cw.exafsdg = np.copy(self.exafsd[cnr])
        
        self.cw.mode = 3 #compare mju
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.energy.append(self.dataClasses[y].energy)
            self.cw.mju.append(self.dataClasses[y].mju)  
        self.cw.plot()        
        
        self.cw.exec_()
        
    def comparexanes(self):
        if self.current < 0:
            return
        
        self.cw = CompareWindow()
        #self.cw.exafs = np.copy(self.exafsd[cnr])
        #self.cw.exafsdg = np.copy(self.exafsd[cnr])
        
        self.cw.mode = 4 #compare xanes
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.E0 = self.dataClasses[y].E0
            self.cw.energy.append(self.dataClasses[y].energy)
            idx = argmin(abs(self.dataClasses[y].energy - self.dataClasses[y].normalizationEnergy))
            
#            self.cw.mju.append(self.mju_bc[y] / self.mju_bc[y][idx]) 
            self.cw.mju.append(self.dataClasses[y].mjuMinusVictoreen / self.dataClasses[y].mju0[idx]) 
        self.cw.plot()        
        
        self.cw.exec_()
        
    def compareexafs(self):
        
        self.cw = CompareWindow()
        #self.cw.exafs = np.copy(self.exafsd[cnr])
        #self.cw.exafsdg = np.copy(self.exafsd[cnr])
        
        self.cw.mode = 0 #compare exafs
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.k.append(self.dataClasses[y].k)
            self.cw.exafs.append(self.dataClasses[y].exafsZLC)  
        self.cw.plot()        
        
        self.cw.exec_()
    
    def compareft(self):
        
        self.cw = CompareWindow()
        #self.cw.exafs = np.copy(self.exafsd[cnr])
        #self.cw.exafsdg = np.copy(self.exafsd[cnr])
        
        self.cw.mode = 1 #compare exafs
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.r.append(self.dataClasses[y].rZLC)
            self.cw.fr.append(self.dataClasses[y].efrZLC)  
            self.cw.fi.append(self.dataClasses[y].efiZLC)
        self.cw.plot()        
        
        self.cw.exec_()
        
    def comparebft(self):
        
        self.cw = CompareWindow()
        #self.cw.exafs = np.copy(self.exafsd[cnr])
        #self.cw.exafsdg = np.copy(self.exafsd[cnr])
        
        self.cw.mode = 2 #compare exafs
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.bftk.append(self.dataClasses[y].bftk)
            self.cw.bftexafs.append(self.dataClasses[y].bftEXAFS)  
        self.cw.plot()        
        
        self.cw.exec_()
        

    
    def compareXes(self):
        self.cw = CompareWindow()
        
        self.cw.mode = 10 #compare XES original
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.energy.append(self.dataClasses[y].energy)
            self.cw.xes.append(self.dataClasses[y].xes)  
        self.cw.plot()        
        
        self.cw.exec_()  
    
    def compareXesAnorm(self):
        self.cw = CompareWindow()
        
        self.cw.mode = 10 #compare XES original
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        selected_items = self.lstSpectra.selectedItems()
        
        if len(selected_indexes) == 0:
            return

        for i in range(len(list(selected_items))):
            self.cw.labels.append(str(selected_items[i].text()))
        
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            self.cw.energy.append(self.dataClasses[y].energy)
            self.cw.xes.append(self.dataClasses[y].xesAreaNorm)  
        self.cw.plot()        
        
        self.cw.exec_()
        
    def onpick3(self, event):
        
        ind = event.ind
        print(ind)
        cnr = self.current
        self.dataClasses[cnr].energy = delete(self.dataClasses[cnr].energy, ind[0])
        self.dataClasses[cnr].mju = delete(self.dataClasses[cnr].mju, ind[0])
        self.dataClasses[cnr].processExpData()
#        self.ees_step1()
#        self.lstSpectraItemClicked()
        
    def redraw_for_remove(self, state):
        cnr = self.current
        self.ax_abs.clear()
        self.ax_abs2.clear()
        self.ax_abs.scatter(self.energy[cnr], self.mju[cnr], s=1, picker=state)
 
        line2, line3 = self.ax_abs.plot(self.energy[cnr], self.mjub[cnr],
                                      self.energy[cnr], self.mju0_mjub[cnr]) 

        line2.set_color('r')
        line3.set_color('r')
        
    def averageMju(self): 

        selected_indexes = self.lstSpectra.selectedIndexes()
        
        if len(selected_indexes) == 0:
            return
        
        selectedRows = [x.row() for x in selected_indexes]
        
        #check if all datasets with the same length
        
        elementsCount = asarray([len(self.dataClasses[x].mju) for x in selectedRows])
        if sum( elementsCount - elementsCount[0] ) != 0: #not all arrays are the same length
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Not all arrays are the same length. Can't average.")
#            msgBox.setInformativeText("Do you want to save your changes?")
            msgBox.setStandardButtons(QtGui.QMessageBox.Yes)
            msgBox.exec_()
            return

        
        average = array(self.dataClasses[selectedRows[0]].mju)
        name = "mju_average " + str(selectedRows[0])
        for i in selectedRows[1:]:
            name = name + " + " + str(i)
            average = average + self.dataClasses[i].mju
        average = average / len(list(selected_indexes))
        
        
        self.lstSpectra.addItem(name)
        
        self.dataClasses.append(xaesa_exafs_class(2)) # mju
        
        self.dataClasses[-1].energy =  self.dataClasses[selectedRows[0]].energy
        self.dataClasses[-1].mju =  average
        
        self.dataClasses[-1].E0 =  self.dataClasses[selectedRows[0]].E0
        self.dataClasses[-1].E1 =  self.dataClasses[selectedRows[0]].E0
        self.dataClasses[-1].E2 =  self.dataClasses[selectedRows[0]].E0
        self.dataClasses[-1].E3 =  self.dataClasses[selectedRows[0]].E0
        self.dataClasses[-1].zeroLineCorr =  self.dataClasses[selectedRows[0]].zeroLineCorr
        
        self.dataClasses[-1].kMin =  self.dataClasses[selectedRows[0]].kMin
        self.dataClasses[-1].kMax =  self.dataClasses[selectedRows[0]].kMax
        self.dataClasses[-1].dk =  self.dataClasses[selectedRows[0]].dk
        self.dataClasses[-1].rMin =  self.dataClasses[selectedRows[0]].rMin
        self.dataClasses[-1].rMax =  self.dataClasses[selectedRows[0]].rMax
        self.dataClasses[-1].dr =  self.dataClasses[selectedRows[0]].dr
        
        self.dataClasses[-1].kPower =  self.dataClasses[selectedRows[0]].kPower
        
        self.dataClasses[-1].rMinBft =  self.dataClasses[selectedRows[0]].rMinBft
        self.dataClasses[-1].rMaxBft =  self.dataClasses[selectedRows[0]].rMaxBft
        self.dataClasses[-1].bftWindowParam =  self.dataClasses[selectedRows[0]].bftWindowParam
        
        self.dataClasses[-1].processExpData()
                            
    def averageExafs(self):
        pass
    
    def averageXes(self):
        selected_indexes = self.lstSpectra.selectedIndexes()
        
        if len(selected_indexes) == 0:
            return
        
        selectedRows = [x.row() for x in selected_indexes]
        
        #check if all datasets with the same length
        
        elementsCount = asarray([len(self.dataClasses[x].xes) for x in selectedRows])
        if sum( elementsCount - elementsCount[0] ) != 0: #not all arrays are the same length
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Not all arrays are the same length. Can't average.")
#            msgBox.setInformativeText("Do you want to save your changes?")
            msgBox.setStandardButtons(QtGui.QMessageBox.Yes)
            msgBox.exec_()
            return
        average = array(self.dataClasses[selectedRows[0]].xes)
        name = self.lstSpectra.item(selectedRows[0]).text() + '_average_of_' + str(len(selectedRows))
        for i in selectedRows[1:]:
            average = average + self.dataClasses[i].xes
        average = average / len(list(selected_indexes))
        
        
        self.lstSpectra.addItem(name)
            
        self.dataClasses.append(xaesa_xes_class())
        self.dataClasses[-1].name = name
        
        self.dataClasses[-1].energy = self.dataClasses[selectedRows[0]].energy
        self.dataClasses[-1].xes = average


        self.dataClasses[-1].E0 = self.dataClasses[selectedRows[0]].E0
        self.dataClasses[-1].E1 = self.dataClasses[selectedRows[0]].E1
        self.dataClasses[-1].E2 = self.dataClasses[selectedRows[0]].E2
        self.dataClasses[-1].E3 = self.dataClasses[selectedRows[0]].E3
        self.dataClasses[-1].eAreaNormMin = self.dataClasses[selectedRows[0]].eAreaNormMin
        self.dataClasses[-1].eAreaNormMax = self.dataClasses[selectedRows[0]].eAreaNormMax
        
        #process opened file
        self.dataClasses[-1].removeBackground()
        self.dataClasses[-1].areaNormalize()
        
        sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
        sel_item.setForeground(QtCore.Qt.darkMagenta)     
        
        
    def kPowerChange(self):
        cnr = self.current
        newKPower = float(self.edtkpow.text())
        self.dataClasses[cnr].changeKPower(newKPower)
        self.lstSpectraItemClicked() #update graphs

    
    def extractParamsChange(self):
        
        startTime = timer()
        
        cnr = self.current
        
        self.dataClasses[cnr].E0 = float(self.edtE0.text())
        self.dataClasses[cnr].E1 = float(self.edtE1.text())
        self.dataClasses[cnr].E2 = float(self.edtE2.text())
        self.dataClasses[cnr].E3 = float(self.edtE3.text())
        self.dataClasses[cnr].kPower = float(self.edtkpow.text())
        self.dataClasses[cnr].zeroLineCorr = float(self.edtsm.text())
        
        self.dataClasses[cnr].normalizationMode = int(self.chkExafsNormalization.isChecked())
        self.dataClasses[cnr].normalizationEnergy = float(self.edtExafsNormEnergy.text())

        self.dataClasses[cnr].kMin = float(self.edtkmin.text())
        self.dataClasses[cnr].kMax = float(self.edtkmax.text())
        self.dataClasses[cnr].dk = float(self.edtdk.text())
        self.dataClasses[cnr].rMin = float(self.edtrmin.text())
        self.dataClasses[cnr].rMax = float(self.edtrmax.text())
        self.dataClasses[cnr].dr = float(self.edtdr.text())
        
        self.dataClasses[cnr].mju0PolinomialDegree = int(self.edtmju0poldegree.text())
        
        self.dataClasses[cnr].rMinBft = float(self.edtrminbft.text())
        self.dataClasses[cnr].rMaxBft = float(self.edtrmaxbft.text())
        self.dataClasses[cnr].bftWindowParam = float(self.edtbftwindowparam.text())
        
        self.dataClasses[cnr].redoExtraction()
        
        endTime = timer()
        self.lblStatus.setText("Extraction done in {:.4f} seconds".format(endTime - startTime) )
            
        self.lstSpectraItemClicked() #update graphs

    
    def ftBftParamsChange(self):
        
        startTime = timer()
        
        cnr = self.current
        
        self.dataClasses[cnr].E0 = float(self.edtE0.text())
        self.dataClasses[cnr].E1 = float(self.edtE1.text())
        self.dataClasses[cnr].E2 = float(self.edtE2.text())
        self.dataClasses[cnr].E3 = float(self.edtE3.text())
        self.dataClasses[cnr].kPower = float(self.edtkpow.text())
        self.dataClasses[cnr].zeroLineCorr = float(self.edtsm.text())
        
        self.dataClasses[cnr].normalizationMode = int(self.chkExafsNormalization.isChecked())
        self.dataClasses[cnr].normalizationEnergy = float(self.edtExafsNormEnergy.text())

        self.dataClasses[cnr].kMin = float(self.edtkmin.text())
        self.dataClasses[cnr].kMax = float(self.edtkmax.text())
        self.dataClasses[cnr].dk = float(self.edtdk.text())
        self.dataClasses[cnr].rMin = float(self.edtrmin.text())
        self.dataClasses[cnr].rMax = float(self.edtrmax.text())
        self.dataClasses[cnr].dr = float(self.edtdr.text())
        
        self.dataClasses[cnr].mju0PolinomialDegree = int(self.edtmju0poldegree.text())
        
        self.dataClasses[cnr].rMinBft = float(self.edtrminbft.text())
        self.dataClasses[cnr].rMaxBft = float(self.edtrmaxbft.text())
        self.dataClasses[cnr].bftWindowParam = float(self.edtbftwindowparam.text())
        
        self.dataClasses[cnr].redoFtBft()
        
        endTime = timer()
        self.lblStatus.setText("FT & BFT done in {:.4f} seconds".format(endTime - startTime) )
            
        self.lstSpectraItemClicked() #update graphs
        
    def bftParamsChange(self):
        
        startTime = timer()
        
        cnr = self.current
        
        self.dataClasses[cnr].E0 = float(self.edtE0.text())
        self.dataClasses[cnr].E1 = float(self.edtE1.text())
        self.dataClasses[cnr].E2 = float(self.edtE2.text())
        self.dataClasses[cnr].E3 = float(self.edtE3.text())
        self.dataClasses[cnr].kPower = float(self.edtkpow.text())
        self.dataClasses[cnr].zeroLineCorr = float(self.edtsm.text())
        
        self.dataClasses[cnr].normalizationMode = int(self.chkExafsNormalization.isChecked())
        self.dataClasses[cnr].normalizationEnergy = float(self.edtExafsNormEnergy.text())

        self.dataClasses[cnr].kMin = float(self.edtkmin.text())
        self.dataClasses[cnr].kMax = float(self.edtkmax.text())
        self.dataClasses[cnr].dk = float(self.edtdk.text())
        self.dataClasses[cnr].rMin = float(self.edtrmin.text())
        self.dataClasses[cnr].rMax = float(self.edtrmax.text())
        self.dataClasses[cnr].dr = float(self.edtdr.text())
        
        self.dataClasses[cnr].mju0PolinomialDegree = int(self.edtmju0poldegree.text())
        
        self.dataClasses[cnr].rMinBft = float(self.edtrminbft.text())
        self.dataClasses[cnr].rMaxBft = float(self.edtrmaxbft.text())
        self.dataClasses[cnr].bftWindowParam = float(self.edtbftwindowparam.text())
        
        self.dataClasses[cnr].redoBft()
        
        endTime = timer()
        self.lblStatus.setText("BFT done in {:.4f} seconds".format(endTime - startTime) )
            
        self.lstSpectraItemClicked() #update graphs

        
    def rebinMju(self):
        if self.current < 0:
            return
        
        de1 = self.xaesaSettings.rebinE1
        de2 = self.xaesaSettings.rebinE1E0
        dk = self.xaesaSettings.rebinE0E3
            
        msgBox = QtGui.QMessageBox()
        msgBox.setText("E1, E0+50, E3 energy values will be used for rebinning. \n\
Step size is defined as:\n        %s eV for E < E1\n        %s eV for E1 < E < E2 \n        dk=%s A-1 for E2 < E < E3 \n\
Proceed ?" %(de1, de2, dk))
#            msgBox.setInformativeText("Do you want to save your changes?")
        msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        result = msgBox.exec_()
        
        if result == QtGui.QMessageBox.No:
            return    

        selected_indexes = self.lstSpectra.selectedIndexes()   

        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            
            self.dataClasses[y].changeToRebinedMjuAveraging(rebinType = 'spline', dE1 = de1, dE2=de2, dK=dk, s = 0)
#            self.dataClasses[y].changeToRebinedMjuRdfSmooth()
            sel_item = self.lstSpectra.item(y)
            sel_item.setText( sel_item.text() + "(rebin)")
            
    def rebinMjuRbfSmooth(self):
        if self.current < 0:
            return
            
        de1 = self.xaesaSettings.rebinE1
        de2 = self.xaesaSettings.rebinE1E0
        dk = self.xaesaSettings.rebinE0E3
            
        msgBox = QtGui.QMessageBox()
        msgBox.setText("E1, E0+50, E3 energy values will be used for rebinning. \n\
Step size is defined as:\n        %s eV for E < E1\n        %s eV for E1 < E < E2 \n        dk=%s A-1 for E2 < E < E3 \n\
Proceed ?" %(de1, de2, dk))
#            msgBox.setInformativeText("Do you want to save your changes?")
        msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        result = msgBox.exec_()
        
        if result == QtGui.QMessageBox.No:
            return    

        selected_indexes = self.lstSpectra.selectedIndexes()
        
        sm = float(self.edtRbfSmooth.text())
        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            
            self.dataClasses[y].changeToRebinedMjuAveraging(rebinType = 'rbf-smooth', dE1 = de1, dE2=de2, dK=dk, s = sm)
#            self.dataClasses[y].changeToRebinedMjuRdfSmooth()
            sel_item = self.lstSpectra.item(y)
            sel_item.setText( sel_item.text() + "(rebin)")
                
        
    def backToOriginalMju(self):
        if self.current < 0:
            return
        
        selected_indexes = self.lstSpectra.selectedIndexes()
        

        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            
            self.dataClasses[y].changeToOriginalMju()
            sel_item = self.lstSpectra.item(y)
            sel_item.setText( sel_item.text().replace('(rebin)', ''))

    def rebinExafs(self):
        if self.current < 0:
            return
            
        msgBox = QtGui.QMessageBox()
        msgBox.setText("E1, E0+50, E3 energy values will be used for rebinning. \n\
Step size is defined as:\n        5.00 eV for E < E1\n        0.25 eV for E1 < E < E2 \n        dk=0.02 A-1 for E2 < E < E3 \n\
Proceed ?")
#            msgBox.setInformativeText("Do you want to save your changes?")
        msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        result = msgBox.exec_()
        
        if result == QtGui.QMessageBox.No:
            return    

        selected_indexes = self.lstSpectra.selectedIndexes()
        

        for i in range(len(list(selected_indexes))):
            y= selected_indexes[i].row()
            
            self.dataClasses[y].changeToRebinedMjuAveraging1()
            sel_item = self.lstSpectra.item(y)
            sel_item.setText( sel_item.text() + "(rebin1)")
            
    def rebinXes(self):
        selected_indexes = self.lstSpectra.selectedIndexes()
        
        if len(selected_indexes) == 0:
            return
        
        selectedRows = [x.row() for x in selected_indexes]
        
        text, ok = QtGui.QInputDialog.getText(self,
                                       "Enter rebin parameters",
                                       "Enter rebin parameters Emin Emax dE",
                                       QtGui.QLineEdit.Normal, "7000 7100 0.01")
    
		
        if ok:
            rebinParams = text.split()
            eMin = float(rebinParams[0])
            eMax = float(rebinParams[1])
            dE = float(rebinParams[2])
        else:
            return
        
        print(rebinParams)
        
        for row in selectedRows:
            energy_scale = self.dataClasses[row].energy
            vals = self.dataClasses[row].xes
            
            #check if all points are increasing
            b = diff(energy_scale) > 0
            
            if not all(b):
                msgBox = QtGui.QMessageBox()
                msgBox.setText("Data points must be increasing, this condition is not satisfied.\n\
Remove decreasing data points automaticly ?")
                msgBox.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
                result = msgBox.exec_()
                if result == QtGui.QMessageBox.No:
                    continue
                if result == QtGui.QMessageBox.Yes: #remove bad data points
                    b = append(b, False)
                    energy_scale = self.dataClasses[row].energy[b]
                    vals = self.dataClasses[row].xes[b]
        
#            new_energy, rebined = rebin_mju(energy_scale, vals, 
#                                         rebinParams[0], rebinParams[1], rebinParams[2])
            
            spl = InterpolatedUnivariateSpline(energy_scale, vals)

            new_energy_scale = arange(eMin, eMax, dE)
            
            new_vals = spl(new_energy_scale)
            
            name = self.lstSpectra.item(row).text() + '_rebin'
            
            self.lstSpectra.addItem(name)
            
            self.dataClasses.append(xaesa_xes_class())
            self.dataClasses[-1].name = name
            
            self.dataClasses[-1].energy = new_energy_scale
            self.dataClasses[-1].xes = new_vals
    
    
            self.dataClasses[-1].E0 = self.dataClasses[row].E0
            self.dataClasses[-1].E1 = self.dataClasses[row].E1
            self.dataClasses[-1].E2 = self.dataClasses[row].E2
            self.dataClasses[-1].E3 = self.dataClasses[row].E3
            self.dataClasses[-1].eAreaNormMin = self.dataClasses[row].eAreaNormMin
            self.dataClasses[-1].eAreaNormMax = self.dataClasses[row].eAreaNormMax
            
            #process opened file
            self.dataClasses[-1].removeBackground()
            self.dataClasses[-1].areaNormalize()
            
            sel_item = self.lstSpectra.item(self.lstSpectra.count()-1)
            sel_item.setForeground(QtCore.Qt.darkMagenta) 

        

    def dwcorrection(self):
        cnr = self.current
        
#        dwfactor = float(askstring("Enter DW factor", "Enter DW factor"))
        num, ok = QtGui.QInputDialog.getDouble(self,
                                               "Enter DW factor",
                                               "Enter DW factor",
                                               value = 0.01, 
                                               max = 1,
                                               min = -1,
                                               decimals = 4)
		
        if ok:
            dwfactor = num
        else:
            return
        
        self.dataClasses[cnr].exafs = self.dataClasses[cnr].exafs * exp(-2*dwfactor*self.dataClasses[cnr].k**2)
#        self.exafsd[cnr] = self.dataClasses[cnr].exafs
        self.dataClasses[cnr].exafsZLC = self.dataClasses[cnr].exafs
        self.dataClasses[cnr].redoFtBft()
        self.lstSpectraItemClicked()
        
    def shiftenergyscale(self):
        
        cnr = self.current
        
#        shift = float(askstring("Enter energy shift", "Enter energy shift"))
        num, ok = QtGui.QInputDialog.getDouble(self,
                                               "Energy shift",
                                               "Enter energy shift",
                                               value = 10, 
                                               max = 10000,
                                               min = -10000,
                                               decimals = 2)
		
        if ok:
            shift = num
        else:
            return
        
        self.energy[cnr] = self.energy[cnr] + shift
        self.lstSpectraItemClicked()
        
    def LinCombination(self):
        if self.current < 0:
            return

        lcw = LCWindow()
        
        lcw.mainform = self
        
        msgBox = QtGui.QMessageBox()
        msgBox.setText('Choose mju or EXAFS?')
        msgBox.addButton(QtGui.QPushButton('Mju'), QtGui.QMessageBox.AcceptRole)
        msgBox.addButton(QtGui.QPushButton('EXAFS'), QtGui.QMessageBox.AcceptRole)
        msgBox.addButton(QtGui.QPushButton('XES'), QtGui.QMessageBox.AcceptRole)
        lcw.mju_exafs = msgBox.exec_()
        
        lcw.initUI()      
        
        lcw.exec_()
        
    def openNexus(self):
        self.fn = self.openaddfiledialog_qtgui()
        if self.fn == []:
            self.lblStatus.setText("No files opened")
            return
        
        self.fn.sort()
        
        fn1 = []
        
        for i in range(0, len(self.fn)):
            head, tail = path.split(str(self.fn[i]))
            fn1.append(tail)
            
        self.currentdir = head
        

        #add all file names to the list
        #and fill arrays with data
        self.lstSpectra.clear()
#        self.lstSpectra.addItems(fn1)
        
        #arrays to store analysed data
        self.init_arrays()
        
        if self.chkMjuOrExafs.isChecked(): #only exafs
            msgBox = QtGui.QMessageBox()
            msgBox.setText('Select k power used?')
            msgBox.addButton(QtGui.QPushButton('0'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('1'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('2'), QtGui.QMessageBox.AcceptRole)
            msgBox.addButton(QtGui.QPushButton('3'), QtGui.QMessageBox.AcceptRole)
            kpower = msgBox.exec_()

        for i in range(0, len(self.fn)):
             self.print_hdf5_file_structure(self.fn[i])           
#            h5file = h5py.File(self.fn[i], "r")
            
            
#            print(list(h5file.keys()))
            
    def print_hdf5_file_structure(self, file_name) :
#    """Prints the HDF5 file structure"""
        file = h5py.File(file_name, 'r') # open read-only
        item = file #["/Configure:0000/Run:0000"]
        self.print_hdf5_item_structure(item)
        file.close()
     
    def print_hdf5_item_structure(self, g, offset='    ') :
#        """Prints the input file/group/dataset (g) name and begin iterations on its content"""
        if   isinstance(g,h5py.File) :
            print(g.file, '(File)', g.name)
     
        elif isinstance(g,h5py.Dataset) :
            print('(Dataset)', g.name, '    len =', g.shape, g.dtype)
     
        elif isinstance(g,h5py.Group) :
            print('(Group)', g.name)
     
        else :
            print('WORNING: UNKNOWN ITEM IN HDF5 FILE', g.name)
            sys.exit ( "EXECUTION IS TERMINATED" )
     
        if isinstance(g, h5py.File) or isinstance(g, h5py.Group) :
#            print(list(g.keys()))
            for key in g.keys() :
                subg = g[key]
                print(offset, key ,"   ", subg.name)
                self.print_hdf5_item_structure(subg, offset + '    ')
                
    def spectradifference(self):
        if self.current < 0:
            return

        selected_indexes = self.lstSpectra.selectedIndexes()

        y1= selected_indexes[0].row()
        y2= selected_indexes[1].row()

        
        self.lstSpectra.addItem("difference"+str(y1) + "-" + str(y2))
        
        self.fill_arrays()
        
        self.energy[-1] = []
        self.mju[-1] = []

        self.k[-1] = self.k[y1]
        self.exafs[-1] = self.exafs[y1] - self.exafs[y2]
        self.exafssm[-1] = self.exafs[y1] - self.exafs[y2]
        self.exafsd[-1] = self.exafs[y1] - self.exafs[y2]
        
        self.E0[-1] = self.E0[y1]
        self.E1[-1] = self.E1[y1]
        self.E2[-1] = self.E2[y1]
        self.E3[-1] = self.E3[y1]
        self.zlinecorr[-1] = self.zlinecorr[y1]

        self.kmin[-1] = self.kmin[y1]
        self.kmax[-1] = self.kmax[y1]
        self.dk[-1] = self.dk[y1]
        self.rmin[-1] = self.rmin[y1]
        self.rmax[-1] = self.rmax[y1]
        self.dr[-1] = self.dr[y1]
        
        self.kpow[-1] = self.kpow[y1]

        self.rminbft[-1] = self.rminbft[y1]
        self.rmaxbft[-1] = self.rmaxbft[y1]
        self.bftwindowparam[-1] = self.bftwindowparam[y1]
        self.mjuorexafs[-1] = 1
        
    def xesExtractRedo(self):
        
        cnr = self.current

        self.current =  self.lstSpectra.currentRow()
        cnr = self.current 
        
        self.dataClasses[cnr].E0 = float(self.edtXesE0.text())
        self.dataClasses[cnr].E1 = float(self.edtXesE1.text())
        self.dataClasses[cnr].E2 = float(self.edtXesE2.text())
        self.dataClasses[cnr].E3 = float(self.edtXesE3.text())
        
        self.dataClasses[cnr].eAreaNormMin = float(self.edtXesEANormMin.text())
        self.dataClasses[cnr].eAreaNormMax = float(self.edtXesEANormMax.text())
        
        self.dataClasses[cnr].removeBackground()
        self.dataClasses[cnr].areaNormalize()
        
        self.lstSpectraItemClicked()
        
    def eventFilter(self, sender, event):
        # this is the function that processes internal drop in notesList
        if event.type() == QtCore.QEvent.ChildRemoved: #drop event
            self.dragStopSelectedRows = [x.row() for x in self.lstSpectra.selectedIndexes()]
            self.dragStartSelectedRows.sort()
            if self.dragStopSelectedRows == self.dragStartSelectedRows: #do nothing
                return 0
#            print("Before sorting", self.dragStartSelectedRows, self.dragStopSelectedRows)
            tmpList = [self.dataClasses[i] for i in self.dragStartSelectedRows]
            self.dragStartSelectedRows.sort(reverse=True)
#            print("after sorting", self.dragStartSelectedRows, self.dragStopSelectedRows, tmpList)
            for i in self.dragStartSelectedRows: #remove elements
                self.dataClasses.pop(i)
            
            self.dragStopSelectedRows.sort()
            self.dataClasses[self.dragStopSelectedRows[0]:self.dragStopSelectedRows[0]] = tmpList
        if event.type() == QtCore.QEvent.ChildAdded: #drag start event
            self.dragStartSelectedRows = [x.row() for x in self.lstSpectra.selectedIndexes()]

        return False # don't actually interrupt anything
    
    
    def onpick(self, event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        indexes = [i for i, j in enumerate(self.leg_lines) if j == legline]
        print(indexes)
        for i in indexes:
            origline = self.lines[i]
            vis = not origline.get_visible()
            origline.set_visible(vis)
#        if self.mode == 1:
#            origline = self.lined1[legline]
#            origline.set_visible(vis)
        
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
            if vis:
                legline.set_alpha(1.0)
            else:
                legline.set_alpha(0.2)
        self.fig.canvas.draw()
                
                
    def adaptOpeningParams(self):
        print('Mju          ', self.chkDataTypeMju.isChecked())
        print('Xes          ', self.chkDataTypeXes.isChecked())
        print('Exafs        ', self.chkDataTypeExafs.isChecked())
        print('Transient    ', self.chkDataTypeTrans.isChecked())
        
        if self.chkDataTypeMju.isChecked():
            self.lblDataX.setText("Energy")
            self.lblDataY.setText("Mju")
            self.gb2.show()
            if self.chkCalcMju.isChecked():
                self.gb3.show()
                self.lblI0Col.show()
                self.lblI1Col.show()
                self.edtI0Col.show()
                self.edtICol.show()
                if self.chkTransOrFluo.isChecked(): #transmission selected
                    self.lblI1Col.setText("I1 column")
                    self.edtFluoCols.hide()
                    self.edtICol.show()
                else: #Fluorescence selected
                    self.edtFluoCols.show()
                    self.edtICol.hide()
                    self.lblI1Col.setText("Fluo columns")
                    
            else:
                self.gb3.hide()
                self.lblI0Col.hide()
                self.lblI1Col.hide()
                self.edtI0Col.hide()
                self.edtICol.hide()
                self.edtFluoCols.hide()
        
        if self.chkDataTypeXes.isChecked():
            self.lblDataX.setText("Energy")
            self.lblDataY.setText("Emission")
            self.gb2.hide()
            self.gb3.hide()
            self.lblI0Col.hide()
            self.lblI1Col.hide()
            self.edtI0Col.hide()
            self.edtICol.hide()
            self.edtFluoCols.hide()
        
        if self.chkDataTypeExafs.isChecked():
            self.lblDataX.setText("k")
            self.lblDataY.setText("EXAFS")
            self.gb2.hide()
            self.gb3.hide()
            self.lblI0Col.hide()
            self.lblI1Col.hide()
            self.edtI0Col.hide()
            self.edtICol.hide()
            self.edtFluoCols.hide()
        
        if self.chkDataTypeTrans.isChecked():
            self.lblDataX.setText("Energy")
            self.lblDataY.setText("Transient")
            self.gb2.hide()
            self.gb3.hide()
            self.lblI0Col.hide()
            self.lblI1Col.hide()
            self.edtI0Col.hide()
            self.edtICol.hide()
            self.edtFluoCols.hide()
            
    def saveGUISettings(self):
        settings = QtCore.QSettings('dev', GUI_SETTINGS_ID)
#        settings.setValue('emin', self.edtEmin.text())
#        settings.setValue('emax', self.edtEmax.text())
#        settings.setValue('de', self.edtdE.text())
#        settings.setValue('reflection', self.edtCrystalReflection.text())
#        settings.setValue('crystal', str(self.cmbCrystal.currentIndex()))
        

    def restoreGUISettings(self):
        settings = QtCore.QSettings('dev', GUI_SETTINGS_ID)
#        self.edtEmin.setText(str(settings.value('emin')))
#        self.edtEmax.setText(str(settings.value('emax')))
#        self.edtdE.setText(str(settings.value('de')))
#        self.edtCrystalReflection.setText(str(settings.value('reflection')))
#        try:
#            self.cmbCrystal.setCurrentIndex(int(str(settings.value('crystal'))))
#        except:
#            self.cmbCrystal.setCurrentIndex(0)

if __name__ == '__main__':
    app = QtGui.QApplication(argv)
    window = MyWindow()
    exit(app.exec_())