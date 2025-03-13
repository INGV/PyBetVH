#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This file is part of pyBetVH.

"""

import os
import shutil
import string
import sys
import wx
import numpy as np
# import matplotlib as mlib
# mlib.use('WX')
import getmaps as gmaps
import globalfunctions as gf
import alphabeta as cab
import ventlocation as vl
import viztool as vt


# CLASS BetFrame
class BetFrame(wx.Frame):

    dflDir, workDir, localDir = gf.setDirs()
  
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, )
        icn = os.path.join(self.workDir, "doc", "icons", "betvh.png")
        # self.SetIcon(wx.Icon(icn, wx.BITMAP_TYPE_ANY))
    
        # main sizers
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox_top = wx.BoxSizer(wx.HORIZONTAL)
        hbox_bot = wx.BoxSizer(wx.HORIZONTAL)
    
        # panel top left   
        self.panelTopLeft = wx.Panel(self, wx.ID_ANY)
        vboxTopLeft = wx.StaticBoxSizer(wx.StaticBox(self.panelTopLeft,
                                        wx.ID_ANY, "VOLCANO SELECTION"),
                                        orient = wx.VERTICAL)
    
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        txt = ("Select PVHA folder settings: ")
        hbox1.Add(wx.StaticText(self.panelTopLeft, wx.ID_ANY, txt, size=(-1, -1)), 
                  0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.TOP, 10)
    
        self.butLoadPVHA = wx.Button(self.panelTopLeft, wx.ID_ANY, "Load PVHA")
        hbox1.Add(self.butLoadPVHA, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 10)
        self.Bind(wx.EVT_BUTTON, self.loadPVHA, self.butLoadPVHA)
    
        vboxTopLeft.Add(hbox1, 0, wx.EXPAND|wx.ALL, 0)
    
        self.panelTopLeft.SetSizer(vboxTopLeft)
        hbox_top.Add(self.panelTopLeft, 1, wx.EXPAND|wx.ALL, 5)
    
        # panel top right 
        self.panelTopRight = wx.Panel(self, wx.ID_ANY)
        vbox_top_right = wx.StaticBoxSizer(wx.StaticBox(self.panelTopRight, 
                                           wx.ID_ANY, "SELECTED VOLCANO"), 
                                           orient = wx.VERTICAL)
    
        hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        hbox4.Add(wx.StaticText(self.panelTopRight, wx.ID_ANY, "Name:"), 
                                0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        self.volname_text = wx.StaticText(self.panelTopRight,
                                          wx.ID_ANY, "")
        hbox4.Add(self.volname_text, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        vbox_top_right.Add(hbox4)
    
        hbox41 = wx.BoxSizer(wx.HORIZONTAL)
        hbox41.Add(wx.StaticText(self.panelTopRight, wx.ID_ANY, 
                                 "Central coordinates (UTM):"), 0, 
                                 wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        self.volc_center_text = wx.StaticText(self.panelTopRight, 
                                              wx.ID_ANY, " ")
        hbox41.Add(self.volc_center_text, 0, 
                   wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        vbox_top_right.Add(hbox41)
    
        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.b_vent = wx.Button(self.panelTopRight, 
                                wx.ID_ANY, "Show vent location")
        hbox5.Add(self.b_vent, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.ventLocation, self.b_vent)
        # self.b_vent.Disable()
        vbox_top_right.Add(hbox5)
    
        hbox_top.Add(self.panelTopRight, 1, wx.EXPAND|wx.ALL, 5)
        self.panelTopRight.SetSizer(vbox_top_right)
        self.panelTopRight.Disable()
    
        # panel bottom left    
        self.panelBotLeft = wx.Panel(self, wx.ID_ANY)
        vbox_bot_left = wx.StaticBoxSizer(wx.StaticBox(self.panelBotLeft,
                                          wx.ID_ANY, "EVENT TREE SELECTION"),
                                          orient=wx.VERTICAL)
        self.tree = wx.TreeCtrl(self.panelBotLeft, wx.ID_ANY)
        self.root = self.tree.AddRoot("N123: Eruption")
           
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.onSelChanged, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_EXPANDED, self.onItemExpanded, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_COLLAPSED, self.onItemCollapsed, self.tree)
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.onActivated, self.tree)
    
        vbox_bot_left.Add(self.tree, 1, wx.EXPAND|wx.ALL, 5)
        hbox_bot.Add(self.panelBotLeft, 1, wx.EXPAND|wx.ALL, 5)
        self.panelBotLeft.SetSizer(vbox_bot_left)
        self.panelBotLeft.Disable()
    
        # panel bottom right (a vertical sizer containing 2 panels)
        vbox_bot_right = wx.BoxSizer(wx.VERTICAL)
           
        # sub-panel 1
        self.panelBotRight = wx.Panel(self, wx.ID_ANY)
        vbox0 = wx.BoxSizer(wx.VERTICAL)
    
        vbox1 = wx.StaticBoxSizer(wx.StaticBox(self.panelBotRight, 
                                  wx.ID_ANY, "ABSOLUTE PROBABILITY"),
                                  orient=wx.VERTICAL)
        vbox1.Add(wx.StaticText(self.panelBotRight, 
                                wx.ID_ANY, "Selected Path:"), 0,
                                #wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
                                wx.LEFT, 5)
        self.display_path = wx.StaticText(self.panelBotRight, 
                                          wx.ID_ANY, " ", size=(250, 100),
                                          style=wx.ALIGN_LEFT|wx.TE_WORDWRAP)
        self.display_path.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT,
                                             wx.FONTSTYLE_NORMAL, 
                                             wx.FONTWEIGHT_NORMAL))
    
        # vbox1.Add(self.display_path, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
        vbox1.Add(self.display_path, 0, wx.LEFT, 5)
    
        vbox2 = wx.StaticBoxSizer(wx.StaticBox(self.panelBotRight,
                                  wx.ID_ANY, "CONDITIONAL PROBABILITY"),
                                  orient=wx.VERTICAL)
        vbox2.Add(wx.StaticText(self.panelBotRight, 
                                wx.ID_ANY, "Selected Node:"), 0, 
                                #wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
                                wx.LEFT, 5)
        self.display_node = wx.StaticText(self.panelBotRight, 
                                          wx.ID_ANY, " ", size=(250, 30), 
                                          style=wx.ALIGN_LEFT|wx.TE_WORDWRAP)
        self.display_node.SetFont(wx.Font(8, wx.FONTFAMILY_DEFAULT, 
                                             wx.FONTSTYLE_NORMAL, 
                                             wx.FONTWEIGHT_NORMAL))
    
        # vbox2.Add(self.display_node, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
        vbox2.Add(self.display_node, 0, wx.LEFT, 5)
    
        vbox0.Add(vbox1, 0, wx.EXPAND|wx.ALL, 5)
        vbox0.Add(vbox2, 0, wx.EXPAND|wx.ALL, 5)
    
        self.panelBotRight.SetSizer(vbox0)
        self.panelBotRight.Disable()
        vbox_bot_right.Add(self.panelBotRight, 0, wx.EXPAND|wx.ALL, 0)
    
        # sub-panel 2 -- info and help buttons    
        self.panel_vis = wx.Panel(self, wx.ID_ANY)
    
        vbox3 = wx.BoxSizer(wx.VERTICAL)
    
        hboxComp = wx.BoxSizer(wx.HORIZONTAL)
        self.bComp = wx.Button(self.panel_vis, wx.ID_ANY,
                                name="Compute", label="COMPUTE")
        self.Bind(wx.EVT_BUTTON, self.calcProb, self.bComp)
        hboxComp.Add(self.bComp, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
        self.bComp.Disable()
        vbox3.Add(hboxComp, 0, wx.EXPAND|wx.ALL, 0)
    
        hbox7 = wx.BoxSizer(wx.HORIZONTAL)
        # self.b_info = wx.Button(self.panel_vis, wx.ID_ANY, "INFO")
        # self.Bind(wx.EVT_BUTTON, self.infoFrame, self.b_info)
        # self.b_help = wx.Button(self.panel_vis, wx.ID_ANY, "HELP")
        # self.Bind(wx.EVT_BUTTON, self.help, self.b_help)
        # hbox7.Add(self.b_info, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
        # hbox7.Add(self.b_help, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
    
        self.b_close = wx.Button(self.panel_vis, wx.ID_ANY, "CLOSE")
        self.Bind(wx.EVT_BUTTON, self.closeBet, self.b_close)
        hbox7.Add(self.b_close, 0, wx.ALIGN_BOTTOM|wx.TOP|wx.LEFT, 5)
    
        vbox3.Add(hbox7, 0, wx.EXPAND|wx.TOP, 10)
        
        self.panel_vis.SetSizer(vbox3)
        vbox_bot_right.Add(self.panel_vis, 0, wx.EXPAND|wx.ALL, 0)
    
        hbox_bot.Add(vbox_bot_right, 0, wx.EXPAND|wx.ALL, 0)
    
        vbox.Add(hbox_top, 0, wx.EXPAND|wx.ALL, 5)
        vbox.Add(hbox_bot, 1, wx.EXPAND|wx.ALL, 5)
    
        self.sb = self.CreateStatusBar()
        # self.sb.SetStatusText("You are running on a "+ gf.getOS() + " platform")
        self.SetSizer(vbox)
        self.SetBackgroundColour("#eeeeee")
        # n_displays = wx.Display.GetCount()
        # display = wx.Display(0)
        # #display = d.GetGeometry()
        # _, _, screen_h, screen_w = display.GetGeometry()
        # screen_scale = display.GetScaleFactor()
        # print(n_displays, screen_h, screen_w, screen_scale)
        # #self.SetSize((0.5*screen_w/n_displays, 0.5*screen_h))
        # self.SetSize((0.5*screen_w*screen_scale, 0.5*screen_h*screen_scale))
        n_displays = wx.Display.GetCount()
        display = wx.Display(0)
        _, _, screen_h, screen_w = display.GetGeometry()
        screen_scale = display.GetScaleFactor()
        #print(n_displays, screen_h, screen_w, screen_scale)
        self.SetSize((int(0.5*screen_h), int(0.5*screen_w)))
        self.Centre()
  
  
    def closeBet(self, event): 
        """
        """    
        #self.Destroy()
        wx.CallAfter(self.Close)
  
  
    def loadPVHA(self, event):
        """
        """    
        global voldir
  
        voldir = gf.selDir(self, event)
        if (voldir == ""):
            return
        else:
      
            if os.path.exists(voldir):
                self.readCfgFile(voldir)
                self.volname_text.SetLabel(volname)
                self.treeCreation()
            else:
                msg = ("ERROR:\n File " + v1def_file + " does not exist. "
                       "Check if the Volcano selection has been made correctly.")
                gf.showErrorMessage(self, msg, "ERROR")
                return
  
  
    def calcProb(self, event):
        """
        """    
    
        pr123, pr4, pr5, pr6, pa, pc, pa_ave, pc_ave = self.calcPost()
        pframe = vt.pyBetVizTool(None, wx.ID_ANY, "pyBetVH - Plotting Tool", 
                                 selected_path, selected_node, nodes_flag,
                                 nodes, volname, dtau, sample, vcx, vcy, imgPath,
                                 dip45, nvents, nsizes, nouts, nareas[iout],
                                 lon[iout], lat[iout], idarea[iout],
                                 xmap_min, xmap_max, ymap_min, ymap_max,
                                 intensity[iout], nint[iout], outcomes[iout], intunits[iout],
                                 geom, par1, par2, par3, par4, par5,
                                 pr123, pr4, pr5, pr6, pa, pc, pa_ave, pc_ave)
        pframe.Show(True)
        pframe.Centre()
        return True
      
  
    def infoFrame(self, event):
        hw = helpWindow.HelpWindow(None, wx.ID_ANY, "pyBetVH Info", "doc/info.html")
        hw.Show(True)
  
      
    def onActivated(self, event):
        """
        not used
        """
        pass
  
  
    def onItemCollapsed(self, event):
        """
        not used
        """
        pass
  
  
    def onItemExpanded(self, event):
        """
        Collapse unselected tree
        """
        item_ex = event.GetItem()
        sel = self.tree.GetItemText(item_ex)
        if self.tree.GetItemParent(item_ex).IsOk():
            par = self.tree.GetItemParent(item_ex)
            item, cookie = self.tree.GetFirstChild(par)
            while item:
                if (item != item_ex):
                    self.tree.Collapse(item)
                item, cookie = self.tree.GetNextChild(par, cookie)
  
  
    def onSelChanged(self, event):
        """
        Get selected node and full path
        """
        
        #global num_sel_nodes
        global selected_outcome
        global selected_node, selected_path
        global nodes_flag
    
        nodes_flag = [1, 1, 1]
        item_sel = event.GetItem()
    
        if self.tree.GetItemParent(item_sel).IsOk():
            # get selected node
            selected_node = self.tree.GetItemText(item_sel)
      
            # get selected full path
            root = self.tree.GetItemText(self.tree.GetRootItem())
            par = self.tree.GetItemParent(item_sel)
            node = self.tree.GetItemText(par)
            path = [selected_node]
            while (node != root):
                path.append(node)
                par = self.tree.GetItemParent(par)
                node = self.tree.GetItemText(par)
      
            path.append(root)
            path.reverse()
            selected_path = str(path[0])
            for i in range(1, len(path)):
                selected_path = selected_path + "\n" + str(path[i])
      
            # re-defining nodes flags for input.dat file 
            for i in range(1, len(path)):
                if (i == 1):  # node 4 -- location
                    for index,item in enumerate(nodes4):
                        if (item == path[i]):
                            nodes_flag.append(index)
                if (i == 2):  # node 5 -- size
                    for index,item in enumerate(nodes5):
                        if (item == path[i]):
                            nodes_flag.append(index + 1)
                if (i == 3):
                    for index,item in enumerate(nodes6):
                        if (item == path[3]):
                            selected_outcome = index
                            nodes_flag.append(index + 1)
                if (i == 4):
                    nodes_flag.append(1)
                if (i == 5):
                    nodes_flag.append(1)
    
        else:
            selected_node = self.tree.GetItemText(self.tree.GetRootItem())
            selected_path = self.tree.GetItemText(self.tree.GetRootItem())
    
        self.display_node.SetLabel(selected_node)
        self.display_path.SetLabel(selected_path)
  
  
    def ventLocation(self, event):
        """
        It opens a new frame showing the vent location by creating a
        new object  vl = ventLocation.VentLocation(), which is defined
        in modules/ventLocation.py
        """    
           
        if (imgPath == ""):
            msg = ("WARNING\nImage map file "
                   "in" + vol_dir + string.lower(vol_name)  + " does not exist."
                   "You can continue but no background map will be plotted.")
            gf.showWarningMessage(self, msg, "WARNING")
        else:
            if (imgPath[-3:] != "png"):
                msg = ("ERROR\nImage extension format" 
                       " must be .png only.")
                gf.showErrorMessage(self, msg, "ERROR")
                return
        
        pars = [par1, par2, par3, par4, par5]
        mapLimits = [xmap_min, xmap_max, ymap_min, ymap_max]
    
        vloc = vl.VentLocation(self, wx.ID_ANY, "Vent Location", imgPath, 
                                       pars, mapLimits, [vcx, vcy])
        # vloc.SetIcon(wx.Icon(os.path.join(self.workDir, "doc", "icons", 
        #                      "vent_location.png"), wx.BITMAP_TYPE_ANY))
  
  
    def readCfgFile(self, *kargs):
        """
        
        """
    
        global volname, dtau, sample
        global vcx, vcy
        global imgPath
        global geom, par1, par2, par3, par4, par5
        global p1, l1, d1
        global p2, l2, d2
        global p3, l3, d3
        global p4, l4, d4
        global p5, l5, d5
        global p6, l6, pds6, pdt6
        global lon, lat, idarea
        global dip45
        global nvents, nsizes, nouts, nareas, nsizesnouts
        global intensity, nint, outcomes, intunits
        global p78file, pd78file
        global xmap_min, xmap_max, ymap_min, ymap_max
    
        volcano = kargs[0]
        
        # Main Settings
        tmp = gf.readMainSettings(volcano)
        volname = tmp[0]
        self.vdir = os.path.join(self.localDir, 
                                 volname.replace(" ", "_"))
        if not os.path.exists(self.vdir):
            os.makedirs(self.vdir)
    
        vctmp = tmp[1].strip().split(",")
        vcx, vcy = float(vctmp[0]), float(vctmp[1])
        geom = tmp[2]
        if (geom == "Field"):
            pars = tmp[3].strip().split(",")
            par1 = float(pars[0])
            par2 = float(pars[1])
            par3 = int(pars[2])
            par4 = int(pars[3])
            par5 = float(pars[4])
        elif (geom == "Cone"):
            pars = tmp[3].strip().split(",")
            par1 = float(pars[0])
            par2 = float(pars[1])
            par3 = float(pars[2])
            par4 = -9999
            par5 = -9999
        else:
            print("ERROR in vent geometry")
                
        utmZone = tmp[4]
        dtau = tmp[5]
        sample = tmp[6]
        bgimage = tmp[7]
        pars = tmp[8].strip().split(",")
        xmap_min = float(pars[0])
        xmap_max = float(pars[1])
        ymap_min = float(pars[2])
        ymap_max = float(pars[3])
        
        # Node 1
        tmp = gf.readNode123(volcano, "Node 1")
        p1 = [tmp[0], 1.0-tmp[0]]
        l1 = int(tmp[1])
        d1 = [tmp[2], tmp[3]-tmp[2]]
        
        # Node 2
        tmp = gf.readNode123(volcano, "Node 2")
        p2 = [tmp[0], 1.0-tmp[0]]
        l2 = int(tmp[1])
        d2 = [tmp[2], tmp[3]-tmp[2]]
        
        # Node 3
        tmp = gf.readNode123(volcano, "Node 3")
        p3 = [tmp[0], 1.0-tmp[0]]
        l3 = int(tmp[1])
        d3 = [tmp[2], tmp[3]-tmp[2]]
        
        # Node 4
        tmp = gf.readNode4(volcano)
        n4file = os.path.join(volcano, tmp[0])
        A = np.loadtxt(n4file)
        p4 = A[:,0]
        d4 = A[:,1]
        l4 = tmp[1]
        nvents = np.shape(A)[0]
        
        # Node 5
        tmp = gf.readNode5(volcano)
        dip45 = tmp[0]
        nsizes = tmp[1]
        n5file = os.path.join(volcano, tmp[2])
        A = np.loadtxt(n5file)
        p5 = np.zeros((nvents,nsizes))
        d5 = np.zeros((nvents,nsizes))
        l5 = np.zeros(nvents)
        
        if (dip45 is False):
            for i in range(nvents):
                p5[i,:] = A[0:nsizes]
                l5[i] = A[nsizes:nsizes+1]
                d5[i,:] = A[nsizes+1:]
        else:
            p5 = A[:,0:nsizes]
            l5 = A[:,nsizes:nsizes+1]
            d5 = A[:,nsizes+1:]
            l5 = np.reshape(l5, (nvents))
        
        # Node 6
        tmp = gf.readNode6(volcano)
        nouts = tmp[0]
        outcomes = tmp[1].split(",")
        intunits = tmp[2].split(",")
        nareastmp = tmp[3].split(",")
        nareas = [int(x) for x in nareastmp]
        n6file = os.path.join(volcano,tmp[4])
        n6intensities = os.path.join(volcano,tmp[5])
        nint = [0]*nouts
        intensity = [0]*nouts
        lon = [0]*nouts
        lat = [0]*nouts
        # imgPath = [0]*nouts
        idarea = [0]*nouts
        
        p6 = [0]*nouts
        l6 = [0]*nouts
        pds6 = [0]*nouts
        pdt6 = [0]*nouts
    
        # xmap_min = 1E8
        # ymap_min = 1E8
        # xmap_max = 0
        # ymap_max = 0
    
        f6 = open(n6intensities, "r")
        for i in range(nouts):
            tmp1 = f6.readline().strip().split()
            nint[i] = len(tmp1)
            intensity[i] = np.array([float(tmp1[j]) for j in range(nint[i])])
      
            # reading file points/areas
            n6ptsareas = os.path.join(volcano, 
                                   tmp[6] + "_out" + str(i+1).zfill(2) + ".txt") 
            A = np.loadtxt(n6ptsareas)
            lon[i] = A[:,0]
            lat[i] = A[:,1]
            idarea[i] = A[:,2]
      
            A = np.loadtxt(n6file + "_out" + str(i+1).zfill(2) + ".txt")
            p6[i] = A[:,0]
            l6[i] = A[:,1]
            pds6[i] = A[:,2]
            pdt6[i] = A[:,3]
            # d6[i] = [pds6, pdt6-pds6]
      
            # tmplim = np.amin(lon[i])
            # if (tmplim < xmap_min):
            #     xmap_min = tmplim
      
            # tmplim = np.amin(lat[i])
            # if (tmplim < ymap_min):
            #     ymap_min = tmplim
      
            # tmplim = np.amax(lon[i])
            # if (tmplim > xmap_max):
            #     xmap_max = tmplim
      
            # tmplim = np.amax(lat[i])
            # if (tmplim > ymap_max):
            #     ymap_max = tmplim
          
        # print(xmap_min, xmap_max, ymap_min, ymap_max)
          
        # background map automatic download
        if (bgimage == "None"):
            if (gf.verifyInternetConn()):
                imgPath = os.path.join(volcano, "map.png")
                gmaps.get_map(xmap_min, ymap_min, xmap_max, ymap_max,
                                  utmZone, imgPath)
            else:
                imgPath = ""
                msg = ("WARNING\nInternet connection is not available.\n"
                       "You can continue but no background map will be plotted.")
                gf.showWarningMessage(self, msg, "WARNING")
    
        else:    
            imgPath = os.path.join(volcano, bgimage)
            if os.path.exists(imgPath):
                pass
            else:
                msg = ("WARNING\nBackground Image path does not exist.\n"
                       "Please check your pybet.cfg file and Load again.")
                gf.showWarningMessage(self, msg, "WARNING")
    
          
        # Node 78
        tmp = gf.readNode78(volcano)
        p78file = [""]*nouts
        pd78file = [""]*nouts
        for i in range(nouts):
            p78file[i] = os.path.join(volcano, 
                                   tmp[0] + "_out" + str(i+1).zfill(2) + ".txt")
            pd78file[i] = os.path.join(volcano, 
                                   tmp[1] + "_out" + str(i+1).zfill(2) + ".txt")
          
        return
      
  
    def calcPost(self, *kargs):
        """ 
        """ 
        
        global nodes, iout
        global pp16, pp3cond, pp4cond, pp5cond, ppcond
    
        # tree selection
        nodes = len(nodes_flag)
        print("Selected Node: {0}".format(nodes))
    
        # selected outcome
        if (nodes >= 6):
            iout = int(nodes_flag[5]-1)
        else:
            iout = 0 
        
        pp_con_ave = []
    
        # Calculating alpha and beta at Nodes 1-6 
        alpha1 = cab.makeAlpha16(2, p1, l1, d1)
        alpha2 = cab.makeAlpha16(2, p2, l2, d2)
        alpha3 = cab.makeAlpha16(2, p3, l3, d3)
        alpha4 = cab.makeAlpha16(nvents, p4, l4, d4)
        
        alpha5 = [0]*(nvents)
        for i in range(nvents):
            alpha5[i] = cab.makeAlpha16(nsizes, p5[i,:], l5[i], d5[i,:])
        
        alpha6 = [0]*(nsizes)
        for i in range(nsizes):
            alpha6[i] = cab.makeAlpha16(2, [p6[iout][i], 1.0-p6[iout][i]], 
                                            l6[iout][i], [pds6[iout][i], 
                                            pdt6[iout][i]-pds6[iout][i]])
        
        # Sampling and Post probability nodes 1-6
        # Nodes 1,2,3
        tmp1 = np.random.dirichlet(alpha1, sample).transpose()
        tmp2 = np.random.dirichlet(alpha2, sample).transpose()
        tmp3 = np.random.dirichlet(alpha3, sample).transpose()
        p123tmp = tmp1[0]*tmp2[0]*tmp3[0] 
        pp123 = 1-np.exp(-p123tmp*dtau)   
        #pp123 = tmp1[0]*tmp2[0]*tmp3[0] # without time window (old method)
        
        pp123_ave = cab.theoreticalAverage(alpha1)*cab.theoreticalAverage(alpha2)*cab.theoreticalAverage(alpha3)
        pp123_ave = 1-np.exp(-pp123_ave*dtau)
    
        # Node 4
        if (nodes < 4):
            ind4 = range(1)
            pp4 = np.ones((sample, len(ind4)))
            pp4_ave = np.ones((len(ind4)))
        else:
            if (nodes_flag[3]==0):
                ind4 = range(nvents) 
            else:
                ind4 = [int(nodes_flag[3]-1)]
            
            tmp4 = np.random.dirichlet(alpha4, sample)
            pp4 = tmp4[:,ind4]
            tmp4_ave = cab.theoreticalAverage(alpha4)
            pp4_ave = cab.theoreticalAverage(alpha4)[ind4]
          
        
        # Node 5
        if (nodes < 5):
            ind5 = range(1)
            pp5 = np.ones((len(ind4), sample, len(ind5)))
            pp5_ave = np.ones((len(ind4), len(ind5)))
        else:
            tmp5 = np.zeros((nvents, sample, nsizes))
            tmp5_ave = np.zeros((nvents, nsizes))
          
            for i in range(nvents):
                tmp = np.random.dirichlet(alpha5[i], sample)
                tmp5[i,:,:] = tmp
                tmp5_ave[i,:] = cab.theoreticalAverage(alpha5[i])
            
            # print(np.shape(tmp))  
            if (int(nodes_flag[4])%2==0):
                tmp = int(nodes_flag[4])/2 - 1
                ind5 = list(range(int(tmp), nsizes)) 
            else:
                ind5 = [int( (nodes_flag[4]+1)/2 - 1 )]
          
            pp5 = np.zeros((len(ind4), sample, len(ind5)))
            pp5_ave = np.zeros((len(ind4), len(ind5)))
            for i in range(len(ind4)):
                pp5[i,:,:] = np.transpose(tmp5[ind4[i], :, ind5])
                pp5_ave[i,:] = np.transpose(tmp5_ave[ind4[i], ind5])
        
        
        # Node 6
        if (nodes < 6):
            ind6 = range(1)
            pp6 = np.ones((sample, len(ind5)))
            pp6_ave = np.ones((len(ind5)))
        else:
            tmp6 = np.zeros((sample, nsizes))
            tmp6_ave = np.zeros((nsizes))
            for i in range(nsizes):
                tmp = np.random.dirichlet(alpha6[i], sample).transpose()
                tmp6[:,i] = tmp[0]
                tmp6_ave[i] = cab.theoreticalAverage(alpha6[i])
        
            pp6 = tmp6[:,ind5]
            pp6_ave = tmp6_ave[ind5]
        
    
        # Absolute probability nodes 1-6
        ind = 0
        pp16 = np.zeros((sample, len(ind4)*len(ind5)))
        pp16_ave = np.zeros((len(ind4)*len(ind5)))
        for i in range(len(ind4)):
            for j in range(len(ind5)):
                pp16[:,ind] = pp123[:]*pp4[:,i]*pp5[i,:,j]*pp6[:,j]
                pp16_ave[ind] = pp123_ave*pp4_ave[i]*pp5_ave[i,j]*pp6_ave[j]
                ind += 1  
        
        # print(np.shape(pp16))
        pp16sumabs = np.sum(pp16, axis=1)
        
    
        # Conditional probability 1-6
        if (nodes == 3):
            pp3cond = pp123
            # pp_con_ave = pp123_ave
            # pp_abs_ave = pp123_ave
            return pp123, pp4, pp5, pp6, pp16, pp3cond, pp123_ave, pp16_ave[0]
    
        elif (nodes == 4):
            pp4cond = tmp4                    # p cond each vent
            pp4cond_ave = tmp4_ave                   
            # ppcond = np.sum(pp4cond, axis=1)  # p cond selected vent(s)
            # pp_con_ave = pp4_ave
            # pp_abs_ave = pp123_ave*pp_con_ave
            return pp123, pp4, pp5, pp6, pp16, pp4cond, pp4cond_ave, pp16_ave
        
        elif (nodes == 5): 
            pp5cond = np.zeros((sample,nsizes))
            # pp = np.zeros((sample,nsizes))
            pp5cond_ave = np.zeros((nsizes))
            for j in range(nsizes):
                tmp0 = 0
                tmp1 = 0
                tmp2 = 0
                tmp3 = 0
                for i in range(len(ind4)):
                    tmp0 += pp4[:, i]
                    tmp1 += pp4[:, i]*tmp5[ind4[i], :, j]
                    tmp2 += pp4_ave[i]
                    tmp3 += pp4_ave[i]*tmp5_ave[ind4[i], j]
        
                # p cond each size / selected vents
                pp5cond[:, j] = tmp1/tmp0
                pp5cond_ave[j] = tmp3/tmp2
      
            # p cond selected size(s) / sel vents
            ppcond = np.sum(pp5cond[:, ind5], axis=1)
      
            return pp123, pp4, pp5, pp6, pp16, pp5cond, pp5cond_ave, pp16_ave
        
    
        elif (nodes >= 6):
            phies = 0
            phies_ave = 0
            for i in range(len(ind4)):
                tmp = 0
                tmp_ave = 0
                for j in range(len(ind5)):
                    tmp += pp5[i, :, j]
                    tmp_ave += pp5_ave[i, j]
                  
                phies += pp4[:, i]*tmp
                phies_ave += pp4_ave[i]*tmp_ave
              
            if (nodes == 6):
                tmp1 = 0
                tmp1_ave = 0
                for i in range(len(ind4)):
                    tmp0 = 0
                    tmp0_ave = 0
                    for j in range(len(ind5)):
                        tmp0 += pp5[i, :, j]*pp6[:, j]
                        tmp0_ave += pp5_ave[i, j]*pp6_ave[j]
                  
                    tmp1 = tmp1 + pp4[:, i]*tmp0
                    tmp1_ave = tmp1_ave + pp4_ave[i]*tmp0_ave
        
                pp6cond = tmp1/phies  
                pp6cond_ave = tmp1_ave/phies_ave  
                return pp123, pp4, pp5, pp6, pp16, pp6cond, pp6cond_ave, pp16_ave
        
       
            if (nodes == 7):
              
                iout = int(nodes_flag[5])-1
                
                p78tmp = np.loadtxt(p78file[iout])
                ppp = p78tmp[:, :-1]
                lamb = p78tmp[:, -1]
                
                nrows, ns = ppp.shape
                # if (ns != nint[iout]):
                #     print("ERROR N. intensity: ", ns, nint[iout])
                #     sys.exit()
                
                narea = int(nrows/(nvents*nsizes))
                print("N. aree: {0}".format(narea))
                
                # past data 7-8
                pastdata = np.zeros(narea*len(ind4)*len(ind5))
                psuc = np.zeros((narea*len(ind4)*len(ind5), nint[iout]))
                puns = np.zeros((narea*len(ind4)*len(ind5), nint[iout]))
        
                if os.path.isfile(pd78file[iout]):
                    d = np.loadtxt(pd78file[iout])
                    npastdata = int(d.shape[0]/(narea+2))
                    print("N. Past Data Node 78: ", npastdata)
                    
                    for i in range(npastdata):
                        irow = (narea+2)*i
                        vent_tmp = int(d[irow])-1
                        size_tmp = int(d[irow+1])-1
                        cond1 = len(np.nonzero(vent_tmp==np.array(ind4))[0])
                        cond2 = len(np.nonzero(size_tmp==np.array(ind5))[0])
                        
                        if (cond1!=0 and cond2!=0):
                            ind = vent_tmp*len(ind5)+size_tmp-ind5[0]
                            # print(ind, (irow+narea+1)*len(ind4)*len(ind5))
                            for j in range(irow+2, irow+narea+2):
                                iarea = j-irow-2
                                index = ind+iarea*len(ind4)*len(ind5)
                        
                                tmp = np.nonzero(d[j]>=intensity[iout])
                                nsuc = len(tmp[0])
                                # print(j, index, ind, iarea, len(ind4), len(ind5))
                        
                                if (nsuc==0):
                                    puns[index, 0] += 1
                                elif (nsuc==nint[iout]):
                                    psuc[index, tmp[0]] += 1
                                else:
                                    psuc[index, tmp[0]] += 1
                                    puns[index, nsuc] += 1
                      
                    print("Past Data Read!")
                
                
                ind = 0
                #hcabs = np.zeros((narea,4,nint[iout]))
                #hccon = np.zeros((narea,4,nint[iout]))
                hcabs = np.zeros((narea, nint[iout], sample))
                hccon = np.zeros((narea, nint[iout], sample))
        
                for ia in range(narea):
                    hcabs_area = np.zeros((nint[iout], sample))
                    hccon_area = np.zeros((nint[iout], sample))
                    ind1 = 0
                    jc = 0
                    for j in ind4:
                        kc = 0
                        tmp1 = 0
                        for k in ind5:
                            irow = ia*nvents*nsizes+j*nsizes+k
                            tmp = cab.makeAlphaBeta78(ppp[irow], lamb[irow])
                            alpha = tmp[0] + psuc[ind]
                            beta = tmp[1] + puns[ind]
                            ind += 1
                            
                            #print("{0}\n{1}\n{2}\n{3}\n{4}".format(alpha, beta, psuc[ind], puns[ind], ind))
                            sample_dir = np.zeros((ns, sample))
                            for i in range(ns):
                                tmp = np.random.dirichlet((alpha[i], beta[i]), sample).transpose()
                                sample_dir[i, :] = tmp[0]
                            
                            hc_sample = np.cumprod(sample_dir, axis=0)
                            hcabs_area += hc_sample*pp16[:, ind1]
                            
                            ind1 += 1
                            tmp1 += hc_sample*pp5[jc, :, kc]*pp6[:, kc]
                            kc += 1 
                        
                        hccon_area += pp4[:, jc]*tmp1
                        jc += 1 
                    
                    hccon_area = hccon_area/phies 
                    
                    hcabs[ia] = hcabs_area
                    hccon[ia] = hccon_area
                  
                # points association to areas
                npts = np.shape(idarea[iout])[0]
                ida = idarea[iout].astype(int)
                hcabs2 = np.zeros((npts, nint[iout], sample))
                hccon2 = np.zeros((npts, nint[iout], sample))
                for i in range(npts):
                    ind = ida[i]-1
                    hcabs2[i, :, :] = hcabs[ind, :, :]
                    hccon2[i, :, :] = hccon[ind, :, :]
        
                # print(np.mean(hcabs2, axis=1))
                pp_con_ave = 0
                pp_abs_ave = 0
                return pp123, pp4, pp5, pp6, hcabs2, hccon2, pp_con_ave, pp_abs_ave
    
        # print(" ")
        # print("Sum P Abs")
        # print("Mean = {:.5f}".format(np.mean(pp16sumabs, axis=0)))
        # print("Perc10 = {:.5f}".format(np.percentile(pp16sumabs,10,axis=0)))
        # print("Perc50 = {:.5f}".format(np.percentile(pp16sumabs,50,axis=0)))
        # print("Perc90 = {:.5f}".format(np.percentile(pp16sumabs,90,axis=0)))
        
        # print(" ")
        # print("P Cond")
        # print("Mean = {:.5f}".format(np.mean(ppcond, axis=0)))
        # print("Perc10 = {:.5f}".format(np.percentile(ppcond,10,axis=0)))
        # print("Perc50 = {:.5f}".format(np.percentile(ppcond,50,axis=0)))
        # print("Perc90 = {:.5f}".format(np.percentile(ppcond,90,axis=0)))
  
  
      
  
    def volcanoSelection(self, event, *kargs):
        """ 
        Selection of the volcano from the selected volcano folder kargs[0].
        Reading of the corresponding info in v1def.txt file.
        """ 
        
        global voldir
        voldir = kargs[0]
        
        if os.path.exists(voldir):
          self.readCfgFile(voldir)
          self.volname_text.SetLabel(volname)
          self.treeCreation()
        else:
          msg = ("ERROR:\n" + voldir + " Directory does not exist. "
                 "Check if the Volcano selection has been made correctly.")
          gf.showErrorMessage(self, msg, "ERROR")
          return
        
  
    def treeCreation(self):  
        """
        """
    
        global selected_node, selected_path
        global nodes4, nodes5, nodes6, nodes78
        global nodes_flag
    
        nodes_flag = [1, 1, 1]
        selected_node = "N123: Eruption"
        selected_path = "N123: Eruption"
        
        # node 4 -- vent locations
        if (nvents):
            nodes4 = ["N4: All locations"]
            for i in range(nvents):
                nodes4.append("N4: Location " + str(i+1))
    
        # node 5 -- size
        if (nsizes):
            nodes5 = []
            for i in range(nsizes):
                if (i != nsizes-1):
                    nodes5.append("N5: Size " + str(i+1))
    
                nodes5.append("N5: Size " + str(i+1) + "+")
    
        else:
            msg = "ERROR:\n Bad size value in v1def.txt"
            gf.showErrorMessage(self, msg, "ERROR") 
            return     
    
        # nodes 6 -- outcomes 
        nodes6 = []
        for i in range(nouts):
            nodes6.append("N6: " + outcomes[i].capitalize())
    
        nodes78 = []
        for i in range(nouts):
            nodes78.append("N78: Hazard curves")
    
    
        self.tree.DeleteAllItems()
        self.root = self.tree.AddRoot("N123: Eruption")
        # creating the GUI tree
        for node4 in nodes4:
            n4 = self.tree.AppendItem(self.root, node4)
            for node5 in nodes5:
                n5 = self.tree.AppendItem(n4, node5)
                for i in range(len(nodes6)):
                    n6 = self.tree.AppendItem(n5, nodes6[i])
                    n7 = self.tree.AppendItem(n6, nodes78[i])
    
        self.panelTopRight.Enable()
        self.panelBotLeft.Enable()
        self.panelBotRight.Enable()
        self.bComp.Enable()
    
        selected_node = self.tree.GetItemText(self.tree.GetRootItem())
        selected_path = self.tree.GetItemText(self.tree.GetRootItem())
    
        self.display_node.SetLabel(selected_node)
        self.display_path.SetLabel(selected_path)
        
        self.Layout()



class pyBetGui(wx.App):
    """
    Instance of the main class
    """
    def OnInit(self):
        frame = BetFrame(None, -1, "PyBetVH")
        frame.Show(True)
        self.SetTopWindow(frame)
        frame.Centre()
        return True


# starting the main gui
if __name__ == "__main__":
    app = pyBetGui(0)
    app.MainLoop()
