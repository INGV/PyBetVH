#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 

This file is part of pyBetVH.

"""

import configparser
import os
import sys
# import urllib2
import urllib.request, urllib.error, urllib.parse
import wx

## global variables
# dflDir = os.path.expanduser("~")
# workDir = os.path.dirname(os.path.realpath(__file__))[:-3]

  

def setDirs():
    """
    """
   
    dflDir = os.path.expanduser("~")
    workDir = os.path.dirname(os.path.realpath(__file__))[:-3]
    localDir = os.path.join(dflDir,".betvh")
    if not os.path.exists(localDir):
        os.makedirs(localDir)
  
    return dflDir, workDir, localDir


def readMainSettings(pvhapath):
    """
    """
    print("Input path: {0}".format(pvhapath))
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath, 'pybet.cfg')
    config.read(filecfg)
    vname = config.get('Main Settings', 'Volcano Name')
    vc = config.get('Main Settings', 'Volcano Center')
    shape = config.get('Main Settings', 'Shape')
    geom = config.get('Main Settings', 'Geometry')
    utm = config.get('Main Settings', 'UTM Zone')
    tw = config.getfloat('Main Settings', 'Time Window')
    sp = config.getint('Main Settings', 'Sampling')
    bg = config.get('Main Settings', 'Background Map')
    bg_lims = config.get('Main Settings', 'Map Limits (m UTM)')
    
    return vname, vc, shape, geom, utm, tw, sp, bg, bg_lims


def readNode123(pvhapath, node):
    """
    """
    
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath, 'pybet.cfg')
    config.read(filecfg)
    p = config.getfloat(node, 'Prior probability')
    l = config.getint(node, 'Equivalent N. Data (Lambda)')
    pds = config.getint(node, 'Past Data (Successes)')
    pdt = config.getint(node, 'Past Data (Total)')
    
    return p, l, pds, pdt
  

def readNode4(pvhapath):
    """
    """
    
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath,'pybet.cfg')
    config.read(filecfg)
    f = config.get('Node 4', 'File Name')
    l = config.getint('Node 4', 'Equivalent N. Data (Lambda)')
    
    return f, l
  

def readNode5(pvhapath):
    """
    """
    
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath,'pybet.cfg')
    config.read(filecfg)
    d45 = config.getboolean('Node 5', 'Node 4-5 Dependence')
    nsizes = config.getint('Node 5', 'N. Sizes')
    f = config.get('Node 5', 'File Name')
    
    return d45, nsizes, f
  

def readNode6(pvhapath):
    """
    """
    
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath,'pybet.cfg')
    config.read(filecfg)
    nout = config.getint('Node 6', 'N. Outcomes')
    outcomes = config.get('Node 6', 'Outcomes')
    units = config.get('Node 6', 'Units')
    na = config.get('Node 6', 'N. Areas')
    f1 = config.get('Node 6', 'File Name')
    f2 = config.get('Node 6', 'File Intensities')
    f3 = config.get('Node 6', 'File Points-Areas')
    
    return nout, outcomes, units, na, f1, f2, f3
  

def readNode78(pvhapath):
    """
    """
    
    config = configparser.RawConfigParser()
    filecfg = os.path.join(pvhapath,'pybet.cfg')
    config.read(filecfg)
    f1 = config.get('Node 78', 'File Name Prior')
    f2 = config.get('Node 78', 'File Name Past Data')
    
    return f1, f2


def selDir(self, event):
    """
    Open a dialog to select a directory path
    """
    dfl_dir = os.path.expanduser("~")
    dlg = wx.DirDialog(self, "Select a directory:", defaultPath=dfl_dir,
                       style=wx.DD_DEFAULT_STYLE
                       #| wx.DD_DIR_MUST_EXIST
                       #| wx.DD_CHANGE_DIR
                       )
    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
    else:
        msg = "WARNING\nYou have NOT selected any directory"
        showWarningMessage(self, msg, "WARNING")
        path = ""
  
    dlg.Destroy()
    return path


def selFile(self, event, *kargs):
    """
    upload_file
    It opens a file dialog, it opens the selected file and
    it returns the corresponding path  
    """
    
    dfl_dir = os.path.expanduser("~")
    dlg = wx.FileDialog(self, message="Upload File", defaultDir=dfl_dir, 
                        defaultFile="", wildcard="*.*", 
                        style=wx.FD_OPEN|wx.FD_CHANGE_DIR)
    
    if (dlg.ShowModal() == wx.ID_OK):
        path = dlg.GetPath()
    else:
        msg = "WARNING\nYou have NOT selected any file"
        showWarningMessage(self, msg, "WARNING")
        path = ""
      
    dlg.Destroy()
    return path
  

def verifyInternetConn():
    try:
        #response = urllib2.urlopen('http://maps.google.com/maps', timeout=3)
        response = urllib.request.urlopen("http://maps.google.com/maps", timeout=3)
        return True
    #except urllib2.URLError as err: pass
    except urllib.error.URLError as err: pass
    return False 


def showWarningMessage(self, *kargs):
    """
    It opens a pop-up dialog showing a warning message.
    """
    dlg = wx.MessageDialog(self, kargs[0], kargs[1], wx.OK|wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()


def showErrorMessage(self, *kargs):
    """
    It opens a pop-up dialog showing an error message.
    """
    dlg = wx.MessageDialog(self, kargs[0], kargs[1], wx.OK|wx.ICON_ERROR)
    dlg.ShowModal()
    dlg.Destroy()


def checkMods():
    """
    """
  
    wxMissingMsg = '''
      ERROR! wxPython module is not installed!
      
      Installation:
      
      Linux
      Most of distributions have their pre-compiled package, so you can 
      easily install it from your preferred package manager. The package 
      name is usually python-wxgtk2.8 (debian-based) or wxpython (redhat-
      based). 
      
      Windows and Mac OSX
      Binaries can be downloaded from http://www.wxpython.org/ in the 
      download area (stable version).
      
      Source code
      In case you need to build it from source code you can follow this 
      http://www.wxpython.org/BUILD-2.8.html
      
      '''
    
    npMissingMsg = '''
      ERROR! numPy module is not installed!
      
      Installation:
      
      Linux
      Most of distributions have their pre-compiled package, so you can 
      easily install it from your preferred package manager. The package 
      name is usually python-numpy or just numpy. 
      
      Windows and Mac OSX
      Binaries can be downloaded from http://sourceforge.net/projects/numpy/ 
      
      '''
    
    mplotMissingMsg = '''
      ERROR! Matplotlib module is not installed!
      
      Installation:
      
      Linux
      Most of distributions have their pre-compiled package, so you can 
      easily install it from your preferred package manager. The package 
      name is usually python-matplotlib or just matplotlib. 
      
      Windows and Mac OSX
      Binaries can be downloaded from 
      http://sourceforge.net/projects/matplotlib/files/matplotlib/ 
      '''
      
    sciMissingMsg = '''
      ERROR! SciPy module is not installed!
      
      Installation:
      
      Linux
      Most of distributions have their pre-compiled package, so you can 
      easily install it from your preferred package manager. The package 
      name is usually python-scipy or just scipy. 
      
      Windows and Mac OSX
      Binaries can be downloaded from 
      http://sourceforge.net/projects/scipy/files/ 
      '''
    
    try:
        import numpy
    except ImportError:
        sys.exit(npMissingMsg)
    
    try:
        import wx
    except ImportError:
        sys.exit(wxMissingMsg)
    
    try:
        import matplotlib
    except ImportError:
        sys.exit(mplotMissingMsg)
    
    try:
        import scipy
    except ImportError:
        sys.exit(sciMissingMsg)
    
    return True
