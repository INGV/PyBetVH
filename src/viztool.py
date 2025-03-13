#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 

This file is part of pyBetVH.

"""

# standard modules
import os
import random
import string
import sys
import globalfunctions as gf
import plotlibs 
import wx
import numpy as np
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.patches import Circle, Wedge, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_wx import _load_bitmap
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib as mpl


# some global plotting settings
mpl.rcParams['xtick.direction']='out'
mpl.rcParams['ytick.direction']='out'
mpl.rcParams['axes.labelsize']='10'
mpl.rcParams['xtick.labelsize']='10'
mpl.rcParams['ytick.labelsize']='10'
mpl.rcParams['legend.fontsize']='10'
mpl.rcParams['axes.titlesize']='10'
mpl.rcParams['font.family']='serif'
mpl.rcParams['font.sans-serif']='Times'


class pyBetVizTool(wx.Frame):
    ''' 
    CLASS pyBetVizTool
  
    '''
  
    # default values
    ptSel = 0            # selected point
    hazSel = 0           # selected hazard model
    staSel = 0           # selected statistic
    probTh = 0.01        # selected probability threshold
    intTh = 1.0          # selected intensity threshold
    tw = 0               # selected time window
    figdpi = 75          # 
    figfmt = "png"       # 
    alpha_pts = 0.75     # 
    
    dflDir, workDir, localDir = gf.setDirs()
    
    # initialization for input.dat file 
    seed = random.uniform(-100000, -10000)   # random seed 
  
    def __init__(self, parent, id, title,
                 sel_path, sel_node, nodes_flag,
                 nodes, volname, dtau, sample, vcx, vcy, imgpath,
                 dip45, nvents, nsizes, nouts, nareas,
                 lon, lat, idarea,
                 xmin_map, xmax_map, ymin_map, ymax_map,
                 iml, nint, outcomes, imt,
                 geom, par1, par2, par3, par4, par5, 
                 p123, p4, p5, p6, pabs, pcon, pcon_ave, pabs_ave):
      
  
        self.sel_path = sel_path
        self.sel_node = sel_node
        self.nodes_flag = nodes_flag
        self.nodes = nodes
        self.volname = volname
        self.dtau = dtau
        self.sample = sample
        self.vcx = vcx/1000
        self.vcy = vcy/1000
        self.imgpath = imgpath
        self.dip45 = dip45
        self.nvents = nvents
        self.nsizes = nsizes
        self.nouts = nouts
        self.nareas = nareas
        self.lon = lon/1000
        self.lat = lat/1000
        self.idarea = idarea
        self.iml = iml
        self.nint = nint
        self.outcomes = outcomes
        self.imt = imt
        self.geom = geom
        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.p123 = p123
        self.p4 = p4
        self.p5 = p5
        self.p6 = p6
        self.pabs = pabs
        self.pcon = pcon
        self.pabs_ave = pabs_ave
        self.pcon_ave = pcon_ave
  
        #print np.shape(pabs)
        #print np.shape(pcon)
        #print np.shape(pabs_ave)
        #print np.shape(pcon_ave)
        #print np.shape(pcon)
        #print np.mean(pcon, axis=0)
        
        #print np.sum(np.mean(pabs,axis=0))
        self.xmin = np.amin(self.lon)
        self.xmax = np.amax(self.lon)
        self.ymin = np.amin(self.lat)
        self.ymax = np.amax(self.lat)
  
        self.xmin_map = xmin_map/1000
        self.xmax_map = xmax_map/1000
        self.ymin_map = ymin_map/1000
        self.ymax_map = ymax_map/1000
  
        self.npts = len(self.idarea)
  
        # self.vt_width = 0.75*screen_w
        # self.vt_height = 0.75*screen_h
  
        dist_vc = np.sqrt( (lon-vcx)**2 + (lat-vcy)**2 )
        self.ptSel = np.argmin(dist_vc)                   # selected point
        self.probTh = int(np.mean(self.pabs*1E8))/1E8     # selected prob th
        self.intTh = int(np.median(self.iml))             # selected int th
  
        title = "Visualization Toolkit"
  
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition)
        icn = os.path.join(self.workDir, "doc", "icons", "plotting_tool.png")
        # self.SetIcon(wx.Icon(icn, wx.BITMAP_TYPE_ANY))
        
        
        if (self.nodes >=3 and self.nodes <= 6):
            self.hbox = wx.BoxSizer(wx.HORIZONTAL)
            self.p1 = wx.Panel(self, wx.ID_ANY)
            self.p1.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, 
                                       wx.FONTWEIGHT_NORMAL))
            
            # left panel
            self.vbox1 = wx.BoxSizer(wx.VERTICAL)
            
            self.vboxAC = wx.StaticBoxSizer(wx.StaticBox(self.p1, wx.ID_ANY, 
                                           "PROBABILITY"), orient=wx.VERTICAL)
            self.hbox_sta = wx.BoxSizer(orient=wx.HORIZONTAL)
            if (self.nodes == 3):
                txtN123 = ("Absolute and Conditional probabability\n"
                           "at Node 123 are equivalent. \n"
                           "")
                self.txtN123 = wx.StaticText(self.p1, wx.ID_ANY, txtN123, size=(-1,-1))
                self.vboxAC.Add(self.txtN123, 0, wx.ALL, 6)
                txt = ("SELECTED PATH/NODE:\n" + sel_path)
            else:
                self.rb1ac = wx.RadioButton(self.p1, wx.ID_ANY, "ABSOLUTE", 
                                            style=wx.RB_GROUP)
                self.rb2ac = wx.RadioButton(self.p1, wx.ID_ANY, "CONDITIONAL")
                self.rb1ac.SetValue(True)
                self.vboxAC.Add(self.rb1ac, 0, wx.TOP, 2)
                self.vboxAC.Add(self.rb2ac, 0, wx.TOP, 2)
                self.Bind(wx.EVT_RADIOBUTTON, self.selAbsConN16, self.rb1ac)
                self.Bind(wx.EVT_RADIOBUTTON, self.selAbsConN16, self.rb2ac)
                txt = ("SELECTED PATH:\n" + sel_path)
  
            self.selection = wx.StaticText(self.p1, wx.ID_ANY, txt, size=(-1,-1))
            self.selection.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, 
                                           wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            # self.vboxAC.Add(self.selection, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.TOP, 10)
            self.vboxAC.Add(self.selection, 0, wx.TOP, 10)
  
            self.vbox1.Add(self.vboxAC, 0, wx.EXPAND|wx.LEFT|wx.TOP, 10)        
      
            self.showTableAbs(self.pabs, self.pabs_ave)
            self.showTableCon(self.pcon, self.pcon_ave)
            self.lcCon.Hide()
  
            # Export table form
            self.expdata = wx.Button(self.p1, wx.ID_ANY, "Save Table", size=(-1,-1))
            self.expdata.SetToolTip(wx.ToolTip("Export "))
            self.Bind(wx.EVT_BUTTON, self.expTab, self.expdata)
            self.vbox1.Add(self.expdata, 0, wx.ALL|wx.EXPAND, 5)
      
            self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
            b_close_pframe = wx.Button(self.p1, wx.ID_ANY, "Close")
            self.Bind(wx.EVT_BUTTON, self.close_pframe, b_close_pframe)
            # self.hbox4.Add(b_close_pframe, 0, wx.TOP|wx.ALIGN_RIGHT, 10)
            self.hbox4.Add(b_close_pframe, 0, wx.TOP, 10)
          
            # self.vbox1.Add(self.hbox4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
            self.vbox1.Add(self.hbox4, 0, wx.ALL, 5)
            self.p1.SetSizer(self.vbox1)
            self.hbox.Add(self.p1, 0, wx.EXPAND|wx.ALL, 5)
  
  
            # right panel
            self.vbox2 = wx.BoxSizer(wx.VERTICAL)
            self.pnlCanvas = wx.Panel(self, wx.ID_ANY)
            self.nb = wx.Notebook(self.pnlCanvas)
            self.pn4 = plotlibs.pn4Canvas(self.nb)
            self.pn5 = plotlibs.pn5Canvas(self.nb)
            self.pn6 = plotlibs.pn6Canvas(self.nb)
            
            if (self.nodes >= 3):
                if (self.nodes >= 3):
                    self.pa = np.sum(self.pabs, axis=1)
                    self.pa_ave = np.sum(self.pabs_ave)
                
                else:
                    self.pa = self.pabs
                    self.pa_ave = self.pabs_ave
                 
  
                self.nb.AddPage(self.pn4, "ECDF")
                self.pn4.plotAbsoluteProb(None, self.pa, self.pa_ave)
  
            if (self.nodes == 3):
                self.pn4bis = plotlibs.pn4Canvas(self.nb)
                self.nb.AddPage(self.pn4bis, "Pie Chart")
                self.pn4bis.plotConditionalProb(None, self.pcon, self.nodes)
            
            if (self.nodes >= 4 and self.nodes_flag[3] == 0):
                self.nb.AddPage(self.pn5, "Vent Map")
                limitsMap = [self.xmin_map, self.xmax_map, self.ymin_map, self.ymax_map]
                limitsFig = [self.xmin, self.xmax, self.ymin, self.ymax]
                pars = [self.par1, self.par2, self.par3, self.par4, self.par5]
                pvents = np.reshape(self.p123,(len(self.p123),1))*self.p4
                self.pn5.showMap(pvents, limitsMap, limitsFig, pars, self.vcx, self.vcy, 
                                 self.imgpath)
                # if (self.dip45 == 2):
                  # self.show_vent_map_abs_button(self.pabs)
                  # if (self.nodes >= 5):
                    # self.nb.AddPage(self.pn6, "Size Map")
                    # self.show_size_map_button(self.pabs)
  
            
            self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onTabChanged)
            # self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.on_tab_changed)
            box_nb = wx.BoxSizer(orient=wx.VERTICAL)
            box_nb.Add(self.nb, 1, wx.EXPAND|wx.ALL, 10)
            self.pnlCanvas.SetSizer(box_nb)
            
            self.vbox2.Add(self.pnlCanvas, 1, wx.EXPAND)
            # self.fig = plt.figure()
            # self.fig.set_dpi(75)
            #self.canvas = FigureCanvas(self.p1, wx.ID_ANY, self.fig)
            #self.toolbar = NavigationToolbar(self.canvas)
            #self.toolbar.DeleteToolByPos(6)
            #self.toolbar.DeleteToolByPos(6)
  
            #self.p2.SetSizer(self.vbox2)
            self.hbox.Add(self.vbox2, 1, wx.EXPAND|wx.ALL, 5)
  
  
        elif (self.nodes == 7):
          
            #hcabs[ia,0] = np.mean(hcabs_area, axis=1) 
            #hccon[ia,0] = np.mean(hccon_area, axis=1) 
            
            # main sizer
            self.hbox = wx.BoxSizer(orient=wx.HORIZONTAL)
            
            # left panel   
            self.pnlLT = wx.Panel(self, wx.ID_ANY)
            self.pnlLT.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, 
                                wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            self.vbox_lt = wx.BoxSizer(orient=wx.VERTICAL)
      
            self.vboxAC = wx.StaticBoxSizer(wx.StaticBox(self.pnlLT, wx.ID_ANY, 
                                           "PROBABILITY"), orient=wx.VERTICAL)
            self.hbox_sta = wx.BoxSizer(orient=wx.HORIZONTAL)
            self.rb1acN78 = wx.RadioButton(self.pnlLT, wx.ID_ANY, "ABSOLUTE", 
                                         style=wx.RB_GROUP)
            self.rb2acN78 = wx.RadioButton(self.pnlLT, wx.ID_ANY, "CONDITIONAL")
            self.rb1acN78.SetValue(True)
            self.vboxAC.Add(self.rb1acN78, 0, wx.TOP, 2)
            self.vboxAC.Add(self.rb2acN78, 0, wx.TOP, 2)
      
            self.Bind(wx.EVT_RADIOBUTTON, self.selAbsConN78, self.rb1acN78)
            self.Bind(wx.EVT_RADIOBUTTON, self.selAbsConN78, self.rb2acN78)
  
            txt = ("SELECTED PATH:\n" + sel_path)
            self.selection = wx.StaticText(self.pnlLT, wx.ID_ANY, txt, size=(-1,-1))
            self.selection.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, 
                                           wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            # self.vboxAC.Add(self.selection, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.TOP, 10)
            self.vboxAC.Add(self.selection, 0, wx.LEFT|wx.TOP, 10)
            self.vbox_lt.Add(self.vboxAC, 0, wx.EXPAND|wx.LEFT|wx.TOP, 10)
      
            self.vboxCP =wx.StaticBoxSizer(wx.StaticBox(self.pnlLT, wx.ID_ANY, 
                                           'CONTROL PANEL'), orient=wx.VERTICAL)
      
            self.vboxCP.Add(wx.StaticText(self.pnlLT, wx.ID_ANY, 
                                 "Area:"), 0, wx.TOP, 6)
            self.hbox_area = wx.BoxSizer(orient=wx.HORIZONTAL)
            self.carea = wx.TextCtrl(self.pnlLT, wx.ID_ANY, size=(80,-1))
            self.carea.SetValue(str(self.ptSel+1))
            # self.hbox_area.Add(self.carea, 0, wx.ALIGN_BOTTOM|wx.TOP, 2)
            self.hbox_area.Add(self.carea, 0, wx.TOP, 2)
            self.barea = wx.Button(self.pnlLT, wx.ID_ANY, 'Update Map', size=(-1,-1))
            self.Bind(wx.EVT_BUTTON, self.selArea, self.barea)
            # self.hbox_area.Add(self.barea, 0, wx.ALIGN_BOTTOM|wx.LEFT, 5)
            self.hbox_area.Add(self.barea, 0, wx.LEFT, 5)
            self.vboxCP.Add(self.hbox_area)
      
            self.vboxCP.Add(wx.StaticText(self.pnlLT, wx.ID_ANY, 
                                 "Probability Threshold (0-1):"), 0, wx.TOP, 6)
            self.hbox_pth = wx.BoxSizer(orient=wx.HORIZONTAL)
            self.cpth = wx.TextCtrl(self.pnlLT, wx.ID_ANY, size=(80,-1))
            self.cpth.SetValue(str(self.probTh))
            self.hbox_pth.Add(self.cpth, 0, wx.TOP, 2)
            self.bpth = wx.Button(self.pnlLT, wx.ID_ANY, 'Update Map', size=(-1,-1))
            self.Bind(wx.EVT_BUTTON, self.selProbabilityTh, self.bpth)
            self.hbox_pth.Add(self.bpth, 0, wx.LEFT, 5)
            self.vboxCP.Add(self.hbox_pth)
      
            intlabel = "Intesity Threshold ({0}-{1}):".format(self.iml[0],self.iml[-1])
            self.vboxCP.Add(wx.StaticText(self.pnlLT, wx.ID_ANY, 
                                 intlabel), 0, wx.TOP, 6)
            self.hbox_ith = wx.BoxSizer(orient=wx.HORIZONTAL)
            self.cith = wx.TextCtrl(self.pnlLT, wx.ID_ANY, size=(80,-1))
            self.cith.SetValue(str(self.intTh))
            self.hbox_ith.Add(self.cith, 0, wx.TOP, 2)
            self.bith = wx.Button(self.pnlLT, wx.ID_ANY, 'Update Map', size=(-1,-1))
            self.Bind(wx.EVT_BUTTON, self.selIntensityTh, self.bith)
            self.hbox_ith.Add(self.bith, 0, wx.LEFT, 5)
            self.vboxCP.Add(self.hbox_ith)
      
            self.stalbl = wx.StaticText(self.pnlLT, wx.ID_ANY, 
                                       "Statistic:")
            self.vboxCP.Add(self.stalbl, 0, wx.TOP, 10)
            self.hbox_sta = wx.BoxSizer(orient=wx.HORIZONTAL)
            self.rb1sta = wx.RadioButton(self.pnlLT, wx.ID_ANY, "Average", 
                                         style=wx.RB_GROUP)
            self.rb2sta = wx.RadioButton(self.pnlLT, wx.ID_ANY, "Percentile")
            self.rb1sta.SetValue(True)
            self.vboxCP.Add(self.rb1sta, 0, wx.TOP, 10)
            self.hbox_sta.Add(self.rb2sta, 0, wx.TOP, 6)
            self.cperc = wx.TextCtrl(self.pnlLT, wx.ID_ANY, size=(40,-1))
            self.cperc.SetValue(str(self.staSel))
            self.hbox_sta.Add(self.cperc, 0, wx.TOP, 2)
            self.bperc = wx.Button(self.pnlLT, wx.ID_ANY, 'Update Map', size=(-1,-1))
            self.Bind(wx.EVT_BUTTON, self.selStatistic, self.bperc)
            self.hbox_sta.Add(self.bperc, 0, wx.LEFT, 5)
            self.vboxCP.Add(self.hbox_sta)
      
            self.vbox_lt.Add(self.vboxCP, 0, wx.EXPAND|wx.ALL, 10) 
            
            
            self.vboxED = wx.StaticBoxSizer(wx.StaticBox(self.pnlLT, wx.ID_ANY, "Export Data"), 
                                        orient=wx.VERTICAL)
      
            # Hazard Map -- Export table form
            #self.tab_hm = wx.StaticBoxSizer(wx.StaticBox(self.pnlLT, wx.ID_ANY, "Hazard Map"), 
                                         #orient=wx.HORIZONTAL)
            #self.fmt_rb1 = wx.RadioButton(self.pnlLT, wx.ID_ANY, ".txt", 
                                          #style=wx.RB_GROUP)
            #self.fmt_rb2 = wx.RadioButton(self.pnlLT, wx.ID_ANY, ".csv")
            #self.fmt_rb3 = wx.RadioButton(self.pnlLT, wx.ID_ANY, ".kml")
            #self.tab_hm.Add(self.fmt_rb1, 0, wx.LEFT|wx.EXPAND, 6)
            #self.tab_hm.Add(self.fmt_rb2, 0, wx.LEFT|wx.EXPAND, 6)
            #self.tab_hm.Add(self.fmt_rb3, 0, wx.LEFT|wx.EXPAND, 6)
            #self.fmt_rb3.Enable(False) 
      
            self.expHM = wx.Button(self.pnlLT, wx.ID_ANY, "Hazard Map", size=(-1,-1))
            self.expHM.SetToolTip(wx.ToolTip("Export "))
            #self.Bind(wx.EVT_BUTTON, lambda event, 
                      #arg=self.zi: self.expHazMapTab(event, arg), self.exptab)
      
            self.Bind(wx.EVT_BUTTON, self.expHazMapTab, self.expHM)
            #self.tab_hm.Add(self.exptab, 0, wx.LEFT|wx.EXPAND, 6)
            self.vboxED.Add(self.expHM, 0, wx.ALL|wx.EXPAND, 5)
      
            # Probability Map -- Export table form
            self.expPM = wx.Button(self.pnlLT, wx.ID_ANY, "Probability Map", size=(-1,-1))
            self.expPM.SetToolTip(wx.ToolTip("Export "))
            self.Bind(wx.EVT_BUTTON, self.expProMapTab, self.expPM)
            self.vboxED.Add(self.expPM, 0, wx.ALL|wx.EXPAND, 5)
      
            # Hazard Curve -- Export table form
            self.expHC = wx.Button(self.pnlLT, wx.ID_ANY, "Hazard Curve", size=(-1,-1))
            self.expHC.SetToolTip(wx.ToolTip("Export "))
            self.Bind(wx.EVT_BUTTON, self.expHazCurTab, self.expHC)
            self.vboxED.Add(self.expHC, 0, wx.ALL|wx.EXPAND, 5)
            
            self.vbox_lt.Add(self.vboxED, 0, wx.EXPAND|wx.ALL, 10) 
  
            self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
            b_close_pframe = wx.Button(self.pnlLT, wx.ID_ANY, "Close")
            self.Bind(wx.EVT_BUTTON, self.close_pframe, b_close_pframe)
            # self.hbox4.Add(b_close_pframe, 0, wx.TOP|wx.ALIGN_RIGHT, 10)
            self.hbox4.Add(b_close_pframe, 0, wx.TOP, 10)
          
            # self.vbox_lt.Add(self.hbox4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
            self.vbox_lt.Add(self.hbox4, 0, wx.ALL, 5)
               
      
            self.pnlLT.SetSizer(self.vbox_lt)
            self.hbox.Add(self.pnlLT, 0, wx.EXPAND|wx.ALL, 5)
            
            # right panel
            vbox_rt = wx.BoxSizer(orient=wx.VERTICAL)
      
            self.pnlCanvas = wx.Panel(self, wx.ID_ANY)
            self.nb = wx.Notebook(self.pnlCanvas)
            self.pn1 = plotlibs.pn1Canvas(self.nb)
            self.pn2 = plotlibs.pn2Canvas(self.nb)
            self.pn3 = plotlibs.pn3Canvas(self.nb)
            
            self.nb.AddPage(self.pn1, "Hazard Map")
            self.nb.AddPage(self.pn2, "Probability Map")
            self.nb.AddPage(self.pn3, "Hazard Curve")
            
            self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onTabChanged)
            #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.on_tab_changed)
            box_nb = wx.BoxSizer(orient=wx.VERTICAL)
            box_nb.Add(self.nb, 1, wx.EXPAND|wx.ALL, 10)
            self.pnlCanvas.SetSizer(box_nb)
            
            vbox_rt.Add(self.pnlCanvas, 1, wx.EXPAND|wx.ALL, 0)
            self.hbox.Add(vbox_rt, 1, wx.EXPAND|wx.ALL, 5)
      
            cid = self.pn1.fig.canvas.mpl_connect('button_press_event', self.onClick)
            self.loadData(self.pabs)
  
        else:
            msg = "ERROR:\n num_sel_nodes must be >= 3 and <= 7"
            gf.showErrorMessage(self, msg, "ERROR") 
            return
      
        self.Bind(wx.EVT_CLOSE, self.onQuit)
        self.SetBackgroundColour("#eeeeee")
        self.sb = self.CreateStatusBar()
        self.sb.SetStatusText("... ")
        self.SetSizer(self.hbox)
        # self.SetSize((self.vt_width, self.vt_height))
        screen_w, screen_h = wx.GetDisplaySize()
        n_displays = wx.Display.GetCount()
        self.SetSize((0.75*screen_w/n_displays, 0.75*screen_h))
        self.Centre()
        #self.Fit()
        self.Show()
        
  
    def savefig(self, event):
        figfmt = self.cfmt.GetValue()
        figdpi = float(self.cdpi.GetValue())
  
        wld = "*." + figfmt
        dlg = wx.FileDialog(self, message="Save Figure as...", 
                            defaultDir=self.dflDir, defaultFile="", 
                            wildcard=wld, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        if (dlg.ShowModal() == wx.ID_OK):
            figpath = dlg.GetPath()
            if (figpath[-4:] == "." + figfmt):
                saved_fig = figpath
            else:
                saved_fig = figpath + "." + figfmt
  
            self.fig.savefig(saved_fig, dpi=figdpi, format=figfmt, 
                           bbox_inches="tight")
          
        dlg.Destroy()
      
  
    def close_pframe(self, event):
        plt.close("all")
        self.Destroy()
  
  
    def expTab(self, event):
        """
        """
  
        if (self.rb1ac.GetValue()):
            values = self.calc_tab_values_abs(self.pabs, 
                                              self.pabs_ave, self.nodes)
        else:
            values = self.calc_tab_values_con(self.pcon, 
                                              self.pcon_ave, self.nodes)
          
        wildcard = "Text files (*.txt; *.dat)|*.txt;*.dat|All files (*.*)|*.*"  
        dlg = wx.FileDialog(self, message="Save File as...", 
                            defaultDir=self.dflDir, defaultFile="", 
                            wildcard=wildcard, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        if (dlg.ShowModal() == wx.ID_OK):
            savepath = dlg.GetPath()
            
            if (savepath[-4:-3] == "."):
                filename = savepath
            else:
                filename = savepath + ".txt"
              
            fp = open(filename, "w")
            if len(values.shape) > 1:
                ncols, nrows = np.shape(values)
                for i in range(nrows):
                    line = " ".join(str(values[j,i]) for j in range(ncols))
                    fp.write("{0}\n".format(line))
  
            else:
                ncols = np.shape(values)[0]
                line = " ".join(str(values[j]) for j in range(ncols))
                fp.write("{0}\n".format(line))
  
                
            fp.close()
        
        dlg.Destroy()
  
  
    def expHazMapTab(self, event):
        """
        """
        wildcard = "Text files (*.txt; *.dat)|*.txt;*.dat| All files (*.*)|*.*"  
        dlg = wx.FileDialog(self, message="Save File as...", 
                            defaultDir=self.dflDir, defaultFile="", 
                            wildcard=wildcard, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        rows = len(self.zi)
        
        if (dlg.ShowModal() == wx.ID_OK):
            savepath = dlg.GetPath()
            
            if (savepath[-4:-3] == "."):
                filename = savepath
            else:
                filename = savepath + ".txt"
              
            fp = open(filename, "w")
            for i in range(rows):
                fp.write("{:f} {:f} {:f}\n".format(self.lon[i], self.lat[i], self.zi[i]))
                
            fp.close()
        
        dlg.Destroy()
      
      
    def expProMapTab(self, event):
        """
        """
          
        wildcard = "Text files (*.txt; *.dat)|*.txt;*.dat|All files (*.*)|*.*"  
        dlg = wx.FileDialog(self, message="Save File as...", 
                            defaultDir=self.dflDir, defaultFile="", 
                            wildcard=wildcard, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        rows = len(self.zp)
        
        if (dlg.ShowModal() == wx.ID_OK):
            savepath = dlg.GetPath()
            
            if (savepath[-4:-3] == "."):
                filename = savepath
            else:
                filename = savepath + ".txt"
              
            fp = open(filename, "w")
            for i in range(rows):
                fp.write("{:f} {:f} {:f}\n".format(self.lon[i], self.lat[i], self.zp[i]))
                
            fp.close()
        
        dlg.Destroy()
  
  
    def expHazCurTab(self, event):
        """
        """
          
        wildcard = "Text files (*.txt; *.dat)|*.txt;*.dat| All files (*.*)|*.*"  
        dlg = wx.FileDialog(self, message="Save File as...", 
                            defaultDir=self.dflDir, defaultFile="", 
                            wildcard=wildcard, style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        rows = len(self.iml)
        
        if (dlg.ShowModal() == wx.ID_OK):
            savepath = dlg.GetPath()
            
            if (savepath[-4:-3] == "."):
                filename = savepath
            else:
                filename = savepath + ".txt"
              
            fp = open(filename, "w")
            for i in range(rows):
                fp.write("{:f} {:f} {:f} {:f} {:f} \n".format(self.iml[i], 
                                                            self.hcp[0][i],
                                                            self.hcp[1][i],
                                                            self.hcp[2][i],
                                                            self.hcp[3][i]))
                
            fp.close()
        
        dlg.Destroy()
  
  
    def loadData(self, *kargs):
        """
        
        """
        self.hc = kargs[0]
        # opening waiting pop-up frame
        # busydlg = wx.BusyInfo("The task is in progress.. please wait", parent=self)
        # wx.Yield()
  
        self.limits = [self.xmin_map, self.xmax_map, self.ymin_map, self.ymax_map]
        
        self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                     self.nareas, self.npts, self.ptSel, 
                                     self.staSel, self.probTh, self.intTh, 
                                     self.imgpath, self.limits, self.iml, 
                                     self.imt)
  
        self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                          self.nareas, self.npts, self.ptSel, 
                                          self.staSel, self.probTh, self.intTh, 
                                          self.imgpath, self.limits, self.iml, 
                                          self.imt)
  
        self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                        self.probTh, self.intTh, self.dtau)
  
  
        #self.sb.SetStatusText("Hazard successfully loaded")
  
  
    def onClick(self, event):
        """
        1) Finding the closest point in the data grid to the point clicked 
           by the mouse on the map at the top canvas. 
        2) Updating the hazard curve plot in the bottom canvas to the 
           selected point.  
        """
        
        if (self.pn1.toolbar.mode == ""):
  
            if (event.inaxes != self.pn1.ax1): 
                return
            else:
                lon1 = min(self.lon)
                lon2 = max(self.lon)
                lat1 = min(self.lat)
                lat2 = max(self.lat)
                xsel, ysel = event.xdata, event.ydata
                if ( xsel >= lon1 and xsel <= lon2 and ysel >= lat1 and ysel <= lat2 ):
                    dist = np.sqrt( (self.lon-xsel)**2 + (self.lat-ysel)**2 )
                    self.ptSel = np.argmin(dist)
                    
                    self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                                 self.nareas, self.npts, self.ptSel, 
                                                 self.staSel, self.probTh, self.intTh, 
                                                 self.imgpath, self.limits, self.iml, 
                                                 self.imt)
              
                    self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                                      self.nareas, self.npts, self.ptSel, 
                                                      self.staSel, self.probTh, self.intTh, 
                                                      self.imgpath, self.limits, self.iml, 
                                                      self.imt)
              
                    self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                                    self.probTh, self.intTh, self.dtau)
  
                    self.carea.SetValue(str(self.ptSel+1))
                  
                else:      
                    return
  
  
    def onQuit(self, event):
        """
        """
        self.Destroy()
  
  
    def onTabChanged(self, event):
        """
        Switching between tabs of bottom canvas
        """
      
        sel = self.nb.GetSelection()
        old = event.GetOldSelection()
        
        #if (sel == 1):
          #self.nb.ChangeSelection(old)
          #msg = ("WARNING:\nThis feature has not been implemented yet")
          #gf.showWarningMessage(self, msg, "WARNING")
          #return
        
        #elif (sel == 2):
          #self.nb.ChangeSelection(old)
          #msg = ("WARNING:\nThis feature has not been implemented yet")
          #gf.showWarningMessage(self, msg, "WARNING")
          #return
        
        #else:
          #pass
  
  
    def selAbsConN16(self, event):
        """
        """
        if (self.rb1ac.GetValue()):
            self.selection.SetLabel("SELECTED PATH:\n" + self.sel_path)
            self.nb.DeleteAllPages()
            self.pn4 = plotlibs.pn4Canvas(self.nb)
            self.pn5 = plotlibs.pn5Canvas(self.nb)
            self.pn6 = plotlibs.pn6Canvas(self.nb)
            self.nb.AddPage(self.pn4, "ECDF")
            self.pn4.plotAbsoluteProb(None, self.pa, self.pa_ave)
            if (self.nodes >= 4 and self.nodes_flag[3] == 0):
                self.nb.AddPage(self.pn5, "Vent Map")
                limitsMap = [self.xmin_map, self.xmax_map, self.ymin_map, self.ymax_map]
                limitsFig = [self.xmin, self.xmax, self.ymin, self.ymax]
                pars = [self.par1, self.par2, self.par3, self.par4, self.par5]
                pvents = np.reshape(self.p123,(len(self.p123),1))*self.p4
                self.pn5.showMap(pvents, limitsMap, limitsFig, pars, self.vcx, self.vcy, 
                               self.imgpath)
                #if (self.dip45 == 2):
                  #self.show_vent_map_abs_button(self.pabs)
                  #if (self.nodes >= 5):
                    #self.nb.AddPage(self.pn6, "Size Map")
                    #self.show_size_map_button(self.pabs)
            self.lcAbs.Show()
            self.lcCon.Hide()
            self.Layout()
                
        else:
            self.selection.SetLabel("SELECTED NODE:\n" + self.sel_node)
            self.nb.DeleteAllPages()
            self.pn4 = plotlibs.pn4Canvas(self.nb)
            self.pn5 = plotlibs.pn5Canvas(self.nb)
            self.pn6 = plotlibs.pn6Canvas(self.nb)
            #nwedges = np.shape(self.pcon)[1]
            
            if (nwedges < 10):
                self.nb.AddPage(self.pn4, "Pie Chart")
                self.pn4.plotConditionalProb(None, self.pcon, self.nodes)
            else:
                pass
                 
            if (self.nodes >= 4 and self.nodes_flag[3] == 0):
                self.nb.AddPage(self.pn5, "Vent Map")
                limitsMap = [self.xmin_map, self.xmax_map, self.ymin_map, self.ymax_map]
                limitsFig = [self.xmin, self.xmax, self.ymin, self.ymax]
                pars = [self.par1, self.par2, self.par3, self.par4, self.par5]
                self.pn5.showMap(self.p4, limitsMap, limitsFig, pars, self.vcx, self.vcy, self.imgpath)
            #if (self.nodes >= 5 and self.nodes_flag[3] == 0):
              #self.show_vent_map_button(self.pcon)
              #if (self.dip45 == 2):
                #self.show_size_map_button(self.pcon)
            #if (self.staSel <= 0 or self.staSel > 100):
              #msg = ("ERROR\nInput value in percentile field is not correct")
              #gf.showErrorMessage(self, msg, "ERROR")
              #self.Raise()
              #return
          
            self.lcCon.Show()
            self.lcAbs.Hide()
            self.Layout()
  
        self.p1.Fit()
          
  
    def selAbsConN78(self, event):
        """
        """
        if (self.rb1acN78.GetValue()):
            self.loadData(self.pabs)
            self.selection.SetLabel("SELECTED PATH:\n" + self.sel_path)
        else:
            self.loadData(self.pcon)
            self.selection.SetLabel("SELECTED NODE:\n" + self.sel_node)
  
        self.pnlLT.Fit()
  
  
    def selArea(self, event):
        """
        """
  
        self.ptSel = int(self.carea.GetValue())
        
        self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                     self.nareas, self.npts, self.ptSel, 
                                     self.staSel, self.probTh, self.intTh, 
                                     self.imgpath, self.limits, self.iml, 
                                     self.imt)
  
        self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                          self.nareas, self.npts, self.ptSel, 
                                          self.staSel, self.probTh, self.intTh, 
                                          self.imgpath, self.limits, self.iml, 
                                          self.imt)
  
        self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                        self.probTh, self.intTh, self.dtau)
  
  
  
  
    def selIntensityTh(self, event):
        """
        """
  
        self.intTh = float(self.cith.GetValue())
  
        self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                     self.nareas, self.npts, self.ptSel, 
                                     self.staSel, self.probTh, self.intTh, 
                                     self.imgpath, self.limits, self.iml, 
                                     self.imt)
  
        self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                          self.nareas, self.npts, self.ptSel, 
                                          self.staSel, self.probTh, self.intTh, 
                                          self.imgpath, self.limits, self.iml, 
                                          self.imt)
  
        self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                        self.probTh, self.intTh, self.dtau)
  
  
  
    def selProbabilityTh(self, event):
        """
        """
  
        self.probTh = float(self.cpth.GetValue())
        self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                     self.nareas, self.npts, self.ptSel, 
                                     self.staSel, self.probTh, self.intTh, 
                                     self.imgpath, self.limits, self.iml, 
                                     self.imt)
  
        self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                          self.nareas, self.npts, self.ptSel, 
                                          self.staSel, self.probTh, self.intTh, 
                                          self.imgpath, self.limits, self.iml, 
                                          self.imt)
  
        self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                        self.probTh, self.intTh, self.dtau)
  
  
    def selStatistic(self, event):
        """
        """
  
        if (self.rb1sta.GetValue()):
            self.staSel = 0
        elif (self.rb2sta.GetValue()):
            self.staSel = int(self.cperc.GetValue())
            if (self.staSel <= 0 or self.staSel > 100):
                msg = ("ERROR\nInput value in percentile field is not correct")
                gf.showErrorMessage(self, msg, "ERROR")
                self.Raise()
                return
            
        else:
            msg = ("ERROR\nInput in Select Statistic is wrong")
            gf.showErrorMessage(self, msg, "ERROR")
            self.Raise()
            return
          
  
        self.zi = self.pn1.hazardMap(self.hc, self.lon, self.lat, self.idarea,
                                     self.nareas, self.npts, self.ptSel, 
                                     self.staSel, self.probTh, self.intTh, 
                                     self.imgpath, self.limits, self.iml, 
                                     self.imt)
  
        self.zp = self.pn2.probabilityMap(self.hc, self.lon, self.lat, self.idarea,
                                          self.nareas, self.npts, self.ptSel, 
                                          self.staSel, self.probTh, self.intTh, 
                                          self.imgpath, self.limits, self.iml, 
                                          self.imt)
  
        self.hcp = self.pn3.hazardCurve(self.hc, self.ptSel, self.iml, self.imt, 
                                        self.probTh, self.intTh, self.dtau)
  
  
  
  ## NOT USED
  #  def export_table(self, event, *args):
  #    s1 = self.fmt_rb1.GetValue()
  #    s2 = self.fmt_rb2.GetValue()
  #    #s3 = str(self.fmt_rb3.GetValue())
  #    
  #    if (s1):
  #      sep = " "
  #      wc = "*.txt"
  #      ext = ".txt"
  #    elif (s2):
  #      sep = ","
  #      wc = "*.csv"
  #      ext = ".csv"
  #    else:
  #      print "error export table"  
  #      
  #    dlg = wx.FileDialog(self, message="Save File as...", 
  #                        defaultDir=self.dflDir, defaultFile="", 
  #                        wildcard=wc, style=wx.SAVE|wx.FD_OVERWRITE_PROMPT)
  #    
  #    rows = len(args[0])
  #    cols = len(args[0][0])
  #    
  #    if (dlg.ShowModal() == wx.ID_OK):
  #      savepath = dlg.GetPath()
  #      
  #      if (savepath[-4:] == ext):
  #        filename = savepath
  #      else:
  #        filename = savepath + ext
  #        
  #      fp = open(filename, "w")
  
  #      for i in range(rows):
  #        for j in range(cols):
  #          if (j == cols-1):
  #            fp.write("%f\n" %(args[0][i][j]))
  #          else:  
  #            fp.write("%f%1s" %(args[0][i][j], sep))
  #          
  #      fp.close()
  #    
  #    dlg.Destroy()
      
  
    #def showTableAbs(self, *kargs):
      #'''
         #Table                
      #'''
      
      #print np.shape(kargs[0])
      #pabs = np.sum(kargs[0], axis=1)
      #values = np.array([ np.percentile(pabs,90,axis=0),
                          #np.percentile(pabs,50,axis=0),
                          #np.percentile(pabs,10,axis=0),
                          #np.mean(pabs, axis=0) ])
  
      ##self.vboxAbs = wx.StaticBoxSizer(wx.StaticBox(self.p1, wx.ID_ANY, 
                                  ##"Absolute Probability"), 
                                  ##orient=wx.VERTICAL)
      #item_list = ("90th Perc", "50th Perc", "10th Perc", "Average")
  
      #self.lcAbs = wx.ListCtrl(self.p1, wx.ID_ANY, 
                       #style=wx.LC_REPORT|wx.LC_NO_HEADER, 
                       #size=(190,120))
      #self.lcAbs.InsertColumn(0, "")
      #self.lcAbs.InsertColumn(1, "")
      #self.lcAbs.SetColumnWidth(0, 100)
      #self.lcAbs.SetColumnWidth(1, 90)
      #context = []
      #tmp_list = []
      #for i in range(len(values)):
        #num_items = self.lcAbs.InsertItem(0, item_list[i])
        #self.lcAbs.SetItem(num_items, 1, "{0:8.4f}".format(values[i]))
  
      ##self.vboxAbs.Add(self.lcAbs, 0, wx.ALL|wx.EXPAND, 5)
      #self.vbox1.Add(self.lcAbs, 0, wx.ALL|wx.EXPAND, 5)
  
  
    def calc_tab_values_abs(self, data, data_ave, n):
        """
        """
        if n == 3:
            values = np.array([data_ave,
                               np.percentile(data, 10),
                               np.percentile(data, 50),
                               np.percentile(data, 90)])
  
        elif n == 4:
            values = np.array([data_ave,
                               np.percentile(data, 10, axis=0),
                               np.percentile(data, 50, axis=0),
                               np.percentile(data, 90, axis=0)])
  
        elif n == 5:
            values = np.array([np.sum(data_ave),
                               np.percentile(np.sum(data, axis=1), 10, axis=0),
                               np.percentile(np.sum(data, axis=1), 50, axis=0),
                               np.percentile(np.sum(data, axis=1), 90, axis=0)])
  
        elif n == 6:
            values = np.array([np.sum(data_ave),
                               np.percentile(np.sum(data, axis=1), 10),
                               np.percentile(np.sum(data, axis=1), 50),
                               np.percentile(np.sum(data, axis=1), 90)])
  
        return values
  
  
    def calc_tab_values_con(self, data, data_ave, n):
        """
        """
        if n == 3 or n == 6:
            values = np.array([data_ave,
                               np.percentile(data, 10),
                               np.percentile(data, 50),
                               np.percentile(data, 90)])
  
        elif n == 4 or n == 5:
            values = np.array([data_ave,
                               np.percentile(data, 10, axis=0),
                               np.percentile(data, 50, axis=0),
                               np.percentile(data, 90, axis=0)])
  
        return values
  
  
  
    def showTableAbs(self, *kargs):
        '''
           Table                
        '''
        
        #pabs = np.sum(kargs[0], axis=1)
        pabs = kargs[0]
        pabs_ave = kargs[1]
        hrow=[]
        
        if (self.nodes == 3):
            header_list = ("", "Eruption")
            hrow.append("Average")
            hrow.append("10th Perc")
            hrow.append("50th Perc")
            hrow.append("90th Perc")
            values = self.calc_tab_values_abs(pabs, pabs_ave, self.nodes)
            ncols = 1
            nrows = len(values)
            nwedges = ncols
  
        elif (self.nodes == 4):
            header_list = ("", "Average", "10th Perc", "50th Perc", "90th Perc")
  #          values = np.array([ pabs_ave,
  #                              np.percentile(pabs,10,axis=0),
  #                              np.percentile(pabs,50,axis=0),
  #                              np.percentile(pabs,90,axis=0) ])
  
            values = self.calc_tab_values_abs(pabs, pabs_ave, self.nodes)
            ncols, nrows = np.shape(values)
            nwedges = nrows
            if (self.nodes_flag[3] == 0):
                for i in range(nrows):
                    hrow.append("Vent " + str(i+1))
            else:
                hrow.append("Vent " + str(self.nodes_flag[3]))
  
        elif (self.nodes == 5):
            
            if (int(self.nodes_flag[4])%2==0):
                tmp = int(self.nodes_flag[4])/2 - 1
                ind5 = range(int(tmp), self.nsizes) 
            else:
                ind5 = [int(self.nodes_flag[4]+1)/2 - 1]
  
            # grouped for sizes
            #if (self.nodes_flag[3] == 0):
              #tmp = np.zeros((self.sample,len(ind5)))
              #for i in range(len(ind5)):
                #tmp[:,i] = np.sum(pabs[:,i:-1:len(ind5)], axis=1) 
            #else:
              #tmp = pabs
  
            #tmp = np.sum(pabs, axis=1)
  
            header_list = ("", "Size")
            hrow.append("Average")
            hrow.append("10th Perc")
            hrow.append("50th Perc")
            hrow.append("90th Perc")
            values = self.calc_tab_values_abs(pabs, pabs_ave, self.nodes)
            ncols = 1
            nrows = len(values)
            nwedges = ncols
  #          nrows = 1
  #          ncols = np.shape(values)[0]
  #          header_list = ("", "Average", "10th Perc", "50th Perc", "90th Perc")
  #          nwedges = nrows
  #          hrow.append("P (Size)")
  #          #for i in ind5:
  #            #hrow.append("Size" + " " + str(i+1))
  
        elif (self.nodes == 6):
            header_list = ("", "Outcome")
            hrow.append("Average")
            hrow.append("10th Perc")
            hrow.append("50th Perc")
            hrow.append("90th Perc")
            values = self.calc_tab_values_abs(pabs, pabs_ave, self.nodes)
            ncols = 1
            nrows = len(values)
            nwedges = ncols
          
  
        else:
            print("Error in showTableAbs")
  
  
        self.lcAbs = wx.ListCtrl(self.p1, wx.ID_ANY, style=wx.LC_REPORT, 
                                 size=(-1,-1))
        context = []
        tmp_list = []
        for i in range(ncols+1):
            self.lcAbs.InsertColumn(i, header_list[i])
        for i in range(nrows):
            self.lcAbs.InsertItem(i, hrow[i])
            for j in range(ncols):
                if (self.nodes == 3 or self.nodes == 5 or self.nodes == 6):
                    self.lcAbs.SetItem(i, j+1, "{0:.3e}".format(values[i]))
                    tmp_list.append(values[j]) 
                else:
                    self.lcAbs.SetItem(i, j+1, "{0:.3e}".format(values[j,i]))
                    tmp_list.append(values[j,i]) 
      
                self.lcAbs.SetColumnWidth(j+1, wx.LIST_AUTOSIZE)
            
            self.lcAbs.SetColumnWidth(0, wx.LIST_AUTOSIZE)
            context.append(tmp_list)
            tmp_list = []

        self.vbox1.Add(self.lcAbs, 0, wx.ALL|wx.EXPAND, 5)
  
  
    def showTableCon(self, *kargs):
        '''
           Table                
        '''
        global nwedges
  
        pcon = kargs[0]
        pcon_ave = kargs[1]
        
        hrow=[]
        if (self.nodes == 3):
            header_list = ("", "Eruption")
            hrow.append("Average")
            hrow.append("10th Perc")
            hrow.append("50th Perc")
            hrow.append("90th Perc")
            values = self.calc_tab_values_con(pcon, pcon_ave, self.nodes)
            ncols = 1
            nrows = len(values)
            nwedges = ncols
  
        elif (self.nodes == 4):
  
            values = self.calc_tab_values_con(pcon, pcon_ave, self.nodes)
            ncols, nrows = np.shape(values)
            header_list = ("", "Average", "10th Perc", "50th Perc", "90th Perc")
            nwedges = nrows
            for i in range(nrows):
                hrow.append("Vent" + " " + str(i + 1))
  
        elif (self.nodes == 5):
  
            values = self.calc_tab_values_con(pcon, pcon_ave, self.nodes)
            ncols, nrows = np.shape(values)
            header_list = ("", "Average", "10th Perc", "50th Perc", "90th Perc")
            nwedges = nrows
            for i in range(nrows):
                hrow.append("Size" + " " + str(i + 1))
  
        elif (self.nodes == 6):
            header_list = ("", "Outcomes")
            hrow.append("Average")
            hrow.append("10th Perc")
            hrow.append("50th Perc")
            hrow.append("90th Perc")
            values = self.calc_tab_values_con(pcon, pcon_ave, self.nodes)
            ncols = 1
            nrows = len(values)
            nwedges = ncols
  
        else:
            print("Error in showTableCon")
  
  
  #      self.vboxCon = wx.StaticBoxSizer(wx.StaticBox(self.p1, wx.ID_ANY, 
  #                                  "Conditional Probability"), 
  #                                  orient=wx.VERTICAL)
  
        self.lcCon = wx.ListCtrl(self.p1, wx.ID_ANY, style=wx.LC_REPORT, size=(-1,-1))
        context = []
        tmp_list = []
        for i in range(ncols+1):
            self.lcCon.InsertColumn(i, header_list[i])
        for i in range(nrows):
            self.lcCon.InsertItem(i, hrow[i])
            for j in range(ncols):
                if (self.nodes == 3 or self.nodes == 6):
                    self.lcCon.SetItem(i, j+1, "{0:.3e}".format(values[i]))
                    tmp_list.append(values[i]) 
                else:
                    self.lcCon.SetItem(i, j+1, "{0:.3e}".format(values[j,i]))
                    tmp_list.append(values[j,i]) 
  
                self.lcCon.SetColumnWidth(j+1, wx.LIST_AUTOSIZE)

            self.lcCon.SetColumnWidth(0, wx.LIST_AUTOSIZE)
            context.append(tmp_list)
            tmp_list = []

        self.vbox1.Add(self.lcCon, 0, wx.ALL|wx.EXPAND, 5)
  
  #      self.vboxCon.Add(self.lcCon, 0, wx.ALL|wx.EXPAND, 5)
  
  #      self.fmt = wx.StaticBoxSizer(wx.StaticBox(self.p1, wx.ID_ANY, 
  #                                   "Format Selection"), 
  #                                   orient=wx.HORIZONTAL)
  #      self.fmt_rb1 = wx.RadioButton(self.p1, wx.ID_ANY, ".txt", 
  #                                    style=wx.RB_GROUP)
  #      self.fmt_rb2 = wx.RadioButton(self.p1, wx.ID_ANY, ".csv")
  #      self.fmt_rb3 = wx.RadioButton(self.p1, wx.ID_ANY, ".kml")
  #      self.fmt.Add(self.fmt_rb1, 0, wx.ALL|wx.EXPAND, 5)
  #      self.fmt.Add(self.fmt_rb2, 0, wx.ALL|wx.EXPAND, 5)
  #      self.fmt.Add(self.fmt_rb3, 0, wx.ALL|wx.EXPAND, 2)
  #      self.fmt_rb3.Enable(False) 
  
  #      self.hbox1.Add(self.fmt, 0, wx.TOP|wx.RIGHT|wx.EXPAND, 5)
  #   
  #      self.exptab1 = wx.Button(self.p1, wx.ID_ANY, "Export", size=(-1,-1))
  #      self.exptab1.SetToolTip(wx.ToolTip("Export Table"))
  #      self.Bind(wx.EVT_BUTTON, lambda event, 
  #                txt=context: self.export_table(event, txt), self.exptab1)
  #      self.hbox1.Add(self.exptab1, 0, wx.TOP|wx.EXPAND, 10)
  
  #      self.sb.Add(self.hbox1, 0, wx.ALL|wx.EXPAND, 5)


