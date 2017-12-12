# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:20:26 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description:
copied & working off of wire.py created by Jorge Chaves
newest changes: changed changed some info included in logs and removed annoying local laptop jazz
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

import os,sys

import wires_constants as const
import wires_functions as func
import wires_toggle as tog

if tog.locallaptop == 0:
    import ROOT
else:
    import wires_local as loc

###############################################################################
# main function
###############################################################################

def wiresMain(ew,lw,subevent,prefix):
    if tog.createlog == 1: func.writeLog_2(timestamp[0], prefix)
    
    allmps_e = []       #list of all early midpoints
    allmps_l = []       #list of all late midpoints
    polyregion_all = [] #list of all polygon regions (as multiple lists of xy coords)
    triareas_all = []   #list of all triangular areas enclosed by wire (of each plane) associated with hit
    badevents = []      #list of subevents that don't meet the badarealimit
    progress = 0        #to show print out of how many tracks finished

    if tog.locallaptop == 0:
        h_angle = ROOT.TH1F("h_angle","",100,-3.15,3.15)
    else:
        h_angle = []
    
    for w_e,w_l,se in zip(ew,lw,subevent):        
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
            if xx[0] != -999:
                if tog.plotEach == 1:
                    plt.plot(xx, yy, const.startcolor, ls='-')
            wires_e.append([xx,yy])
            
            coords = func.getWireCoords(aa, w_l)
            xx = coords[0]
            yy = coords[1]
            if xx[0] != -999:
                if tog.plotEach == 1:
                    plt.plot(xx, yy, const.endcolor, ls='-')
            wires_l.append([xx,yy])
        
        (mp_e, tm_e) = func.getMPorIntersection(wires_e)
        (mp_l, tm_l) = func.getMPorIntersection(wires_l)

        if tm_e == 0:
            ints_e = func.getBigTriangle(wires_e)
        else:
            ints_e = func.getHexagon(mp_e)
        if tm_l == 0:
            ints_l = func.getBigTriangle(wires_l)
        else:
            ints_l = func.getHexagon(mp_l)
        
        ints_all = func.combineAllPoints(ints_e, ints_l)
        ints_all = func.convex_hull(ints_all)
        ints_all = np.array(ints_all)
        
        pa_e = func.PolygonArea(ints_e)
        pa_l = func.PolygonArea(ints_l)
        
        triareas_all.append(pa_e)
        triareas_all.append(pa_l)
        
        if tog.badarealimit != -1 and (func.PolygonArea(ints_e) > tog.badarealimit or func.PolygonArea(ints_l) > tog.badarealimit):
            badevents.append(se)
        else:
            allmps_e.append(mp_e)
            allmps_l.append(mp_l)
            polyregion_all.append(ints_all)
        
        if tog.plotEach == 1:
            if tog.plt_shading == 1:
                ax1 = fig1.add_subplot(111, aspect='equal')
                ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=3*const.hotshading, lw=0))
                ax1.set_xlim(-200,200)
                ax1.set_ylim(-200,200)
                
            if tog.plt_dots == 1: func.plotMPs(mp_e, mp_l)

            if tog.locallaptop == 1:
                plt.savefig(prefix+str(se)+'_'+str(tog.badarealimit)+'.pdf', bbox_inches='tight')
            else:
                plt.savefig('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/events/'+prefix+str(se)+'_'+str(tog.badarealimit)+'.pdf', bbox_inches='tight')
                plt.savefig('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/events/'+prefix+str(se)+'_'+str(tog.badarealimit)+'.png', bbox_inches='tight')
            plt.close()

        if tog.showprogress == 1:
            if tog.plotEach == 1:
                if progress % 25 == 0: print(prefix + str(se) + ".pdf saved\t~ " + str(100*progress/len(subevent)) + "% done")
            else:
                if progress % 25 == 0: print("\t~ " + str(100*progress/len(subevent)) + "% done")
            progress = progress + 1
        
        if tog.createlog == 1: func.writeLog_3(timestamp[0], se, w_e, w_l, mp_e, mp_l, tm_e, tm_l, pa_e, pa_l)
        
    func.plotAllRegions([allmps_e, allmps_l], polyregion_all, prefix)
    
    for gmp_e,gmp_l in zip(allmps_e, allmps_l):
        ang = func.getTrackAngle(gmp_e,gmp_l)
        if ang != -666:
            if tog.locallaptop == 0:
                h_angle.Fill(ang)
            else:
                h_angle.append(ang)

    if tog.locallaptop == 0:
        c1 = ROOT.TCanvas("c1","",900,600)
        h_angle.Draw()
        c1.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c1.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')

    if tog.createlog == 1: func.writeLog_4(timestamp[0], len(subevent), int(np.max(subevent)), np.max(triareas_all), np.average(triareas_all), np.median(triareas_all), len(badevents))
    
    if tog.returnoption == 0:
        return polyregion_all
    elif tog.returnoption == 1:
        return triareas_all
    elif tog.returnoption == 2:
        return badevents
    else:
        return 0

###############################################################################
# data
###############################################################################

# [x,u,v]
ew_beam = []
lw_beam = []
subevent_beam = []

ew_cosmics = []
lw_cosmics = []
subevent_cosmics = []

if tog.locallaptop == 0:
    # Read from ROOT file
    if len(sys.argv) < 2:
        raise ValueError('Missing input file.')
        
    f = ROOT.TFile(sys.argv[1])
    t = ROOT.TTree()
    f.GetObject("tracks",t)
    os.system('mkdir -p plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/events')
    
    ev = 1
    for e in t:
        for i in xrange(0,e.nmatches+e.nmatches2D):
            if e.beam[i] == 1:
                ew_beam.append([int(e.first_hit_X[i]),int(e.first_hit_U[i]),int(e.first_hit_V[i])])
                lw_beam.append([int(e.last_hit_X[i]),int(e.last_hit_U[i]),int(e.last_hit_V[i])])
                subev = '%d.%d' %(e.Event,i)
                subevent_beam.append(float(subev))
            elif e.beam[i] == 0:
                ew_cosmics.append([int(e.first_hit_X[i]),int(e.first_hit_U[i]),int(e.first_hit_V[i])])
                lw_cosmics.append([int(e.last_hit_X[i]),int(e.last_hit_U[i]),int(e.last_hit_V[i])])
                subev = '%d.%d' %(e.Event,i)
                subevent_cosmics.append(float(subev))
            else:
                print("uh oh... something is wrong")
else:
    #tracks_2D_12100.root 
    ew_beam = loc.ew_beam
    lw_beam = loc.lw_beam
    subevent_beam = loc.subevent_beam
    ew_cosmics = loc.ew_cosmics
    lw_cosmics = lw_cosmics
    subevent_cosmics = loc.subevent_cosmics

###############################################################################
# run
###############################################################################

timestamp = func.getTime()

if tog.createlog == 1:
    func.writeLog_1(timestamp, ew_beam, lw_beam, subevent_beam, ew_cosmics, lw_cosmics, subevent_cosmics)

ans1 = wiresMain(ew_beam,lw_beam,subevent_beam,'plot_beam_wires_')
ans2 = wiresMain(ew_cosmics,lw_cosmics,subevent_cosmics,'plot_cosmics_wires_')
