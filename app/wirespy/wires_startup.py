# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:51:20 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: initial setup that gets info & initial plots without any limits
"""

import os,sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import wires_toggle as tog
import wires_functions as func
import wires_constants as const

if tog.locallaptop == 0:
    import ROOT
else:
    import local_wires as loc

from wires_functions import trackInfo, allTrackInfo



###############################################################################
# initial start/setup
###############################################################################

# [x,u,v]
ew_beam = []
lw_beam = []
subevent_beam = []

ew_cosmics = []
lw_cosmics = []
subevent_cosmics = []

n_events = -666

if tog.locallaptop == 0:
    # Read from ROOT file
    if len(sys.argv) < 2:
        raise ValueError('Missing input file.')
        
    f = ROOT.TFile(sys.argv[1])
    t = ROOT.TTree()
    f.GetObject("tracks",t)
    inputfile = sys.argv[1]    
    ev = 1
    n_events = t.GetEntries()
    for e in t:
        for i in xrange(0,e.nmatches+e.nmatches2D):
            if e.beam[i] == 1:
                ew_beam.append([int(e.first_hit_X[i]),int(e.first_hit_U[i]),int(e.first_hit_V[i])])
                lw_beam.append([int(e.last_hit_X[i]),int(e.last_hit_U[i]),int(e.last_hit_V[i])])
                subevent_beam.append([e.Event,i])
            elif e.beam[i] == 0:
                ew_cosmics.append([int(e.first_hit_X[i]),int(e.first_hit_U[i]),int(e.first_hit_V[i])])
                lw_cosmics.append([int(e.last_hit_X[i]),int(e.last_hit_U[i]),int(e.last_hit_V[i])])
                subevent_cosmics.append([e.Event,i])
            else:
                print("uh oh... something is wrong")
else:
    inputfile = "mumbojumbo/test.root"
    ew_beam = loc.ew_beam
    lw_beam = loc.lw_beam
    subevent_beam = loc.subevent_beam
    ew_cosmics = loc.ew_cosmics
    lw_cosmics = loc.lw_cosmics
    subevent_cosmics = loc.subevent_cosmics

inputfile = inputfile.split('/')[-1]

###############################################################################

os.system('mkdir -p plots/'+inputfile.replace('.root',''))
if tog.plotEach == 1:
    os.system('mkdir -p plots/'+inputfile.replace('.root','')+'/events')

timestamp = func.getTime()

func.writeWireID_list(timestamp, inputfile, ew_beam, lw_beam, subevent_beam, ew_cosmics, lw_cosmics, subevent_cosmics)

###############################################################################

def getAllTrackInfo(n_events,ew,lw,subevent,ts,inf,pf):
    alltracks = allTrackInfo(ts, pf, inf, -1, n_events)
    
    func.writeAllTrackInfo_1(alltracks)
    
    if tog.locallaptop == 0 and tog.plt_angles == 1:
        h_angle = ROOT.TH1F("h_angle","",100,-3.15,3.15)
    
    for w_e,w_l,se in zip(ew,lw,subevent):
        newtrack = trackInfo(se=se, w_e=w_e, w_l=w_l)
        wires_e = []    #list of wires associated with early hit
        wires_l = []    #list of wires associated with late hit
        
        if tog.plotEach == 1:
            #plot polygon
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(
                patches.RegularPolygon(
                    (0, 0),      # (x,y)
                    6,           # number of vertices
                    const.radius,      # radius
                    fill = False
                )
            )
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
            if const.plotaxis == 1: func.plotAxis()
        
        #plot 3 wire planes
        for aa in range(3):
            coords = func.getWireCoords(aa, w_e)
            xx = coords[0]
            yy = coords[1]
            if xx[0] != -999 and tog.plotEach == 1:
                plt.plot(xx, yy, const.startcolor, ls='-')
            wires_e.append([xx,yy])
            
            coords = func.getWireCoords(aa, w_l)
            xx = coords[0]
            yy = coords[1]
            if xx[0] != -999 and tog.plotEach == 1:
                plt.plot(xx, yy, const.endcolor, ls='-')
            wires_l.append([xx,yy])
        
        (mp_e, tm_e) = func.getMPorIntersection(wires_e)
        (mp_l, tm_l) = func.getMPorIntersection(wires_l)
        
        newtrack.mp_e = mp_e
        newtrack.mp_l = mp_l
        newtrack.mt_e = sum(tm_e)
        newtrack.mt_l = sum(tm_l)
        
        ang = func.getTrackAngle(mp_e,mp_l)
        newtrack.ang = ang
        if ang != -666 and tog.locallaptop == 0 and tog.plt_angles == 1:
            h_angle.Fill(ang)
        
        if sum(tm_e) == 3:
            ints_e = func.getBigTriangle(wires_e)
        elif sum(tm_e) == 2:
            ints_e = func.getHexagon(mp_e)
        else:
            ints_e = -666
        if sum(tm_l) == 3:
            ints_l = func.getBigTriangle(wires_l)
        elif sum(tm_l) == 2:
            ints_l = func.getHexagon(mp_l)
        else:
            ints_l = -666
        
        ints_e = func.combineAllPoints(ints_e)
        ints_e = func.convex_hull(ints_e)
        ints_l = func.combineAllPoints(ints_l)
        ints_l = func.convex_hull(ints_l)

        newtrack.ic_e = ints_e
        newtrack.ic_l = ints_l
        
        newtrack.ia_e = round(func.PolygonArea(ints_e), const.n_round)
        newtrack.ia_l = round(func.PolygonArea(ints_l), const.n_round)
        
        ints_all = ints_e + ints_l

        ints_all = func.combineAllPoints(ints_all)
        ints_all = func.convex_hull(ints_all)
        
        newtrack.pa = ints_all
        newtrack.tl = round(func.getTrackLength(ints_e, ints_l), const.n_round)
        
        if tog.plotEach == 1:
            #plot individual track info
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=3*const.hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
                
            func.plotMPs(mp_e, mp_l)
            
            dirname = "plots/"+inf.replace('.root','').split('/')[-1]+"/events/"            
            plotname = "plot_"+pf+"_"+str(se[0])+"_"+str(se[1])
            
            if tog.saveplotformat == 0 or tog.saveplotformat == 2:
                plt.savefig(dirname+plotname+'.png', bbox_inches='tight')
                plt.close()
            if tog.saveplotformat == 1 or tog.saveplotformat == 2:
                plt.savefig(dirname+plotname+'.pdf', bbox_inches='tight')
                plt.close()
        
        #if newtrack.ia_e < 0.433 or newtrack.ia_l < 0.433:
        alltracks.add_track(newtrack)
        
    if tog.locallaptop == 0 and tog.plt_angles == 1:
        c1 = ROOT.TCanvas("c1","",900,600)
        h_angle.Draw()
        plotname = "plots/"+inf.replace('.root','').split('/')[-1]+"/h_angle_plot_"+pf
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c1.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c1.Print(plotname+'.png')
    
    func.writeAllTrackInfo_3(alltracks, -1)
    
    func.writeTrackInfo(alltracks)
    if tog.plt_lines == 1 or tog.plt_shading == 1:
        func.plotAllRegions(alltracks)
    
    return alltracks

#

###############################################################################
# run
###############################################################################

tracks_beam = getAllTrackInfo(n_events, ew_beam,lw_beam,subevent_beam,timestamp,inputfile,'beam')
tracks_cosmics = getAllTrackInfo(n_events, ew_cosmics,lw_cosmics,subevent_cosmics,timestamp,inputfile,'cosmics')
