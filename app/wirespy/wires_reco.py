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
import math
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

def wiresMain(n_events, ec,lc,subevent,prefix):
    if tog.createlog == 1: func.writeLog_2(timestamp[0], prefix)
    
    allmps_e = []       #list of all early midpoints
    allmps_l = []       #list of all late midpoints
    polyregion_all = [] #list of all polygon regions (as multiple lists of xy coords)
    triareas_all = []   #list of all triangular areas enclosed by wire (of each plane) associated with hit
    badevents = []      #list of subevents that don't meet the badarealimit
    progress = 0        #to show print out of how many tracks finished
    l_intercept = []

    if tog.locallaptop == 0:
        h_angle = ROOT.TH1F("h_angle","",100,-3.15,3.15)
        h_2D_track_length = ROOT.TH1F("h_2D_track_length","",100,0,300)
        h_3D_track_length = ROOT.TH1F("h_3D_track_length","",100,0,300)
        h_area_density_1 = ROOT.TH1F('h_area_density_1','',6,0,6)
        h_area_density_2 = ROOT.TH1F('h_area_density_2','',2,0,2)
    else:
        h_angle = []
        h_2D_track_length = []
        h_3D_track_length = []
        h_area_density_1 = []
        h_area_density_2 = []
    
    for c_e,c_l,se in zip(ec,lc,subevent):        
        
        x1 = const.radius*math.cos(math.pi/6)/500. * c_e[0]
        y1 = const.radius*math.cos(math.pi/6)/500. * c_e[1]
        z1 = const.radius*math.cos(math.pi/6)/500. * c_e[2]
        x2 = const.radius*math.cos(math.pi/6)/500. * c_l[0]
        y2 = const.radius*math.cos(math.pi/6)/500. * c_l[1]
        z2 = const.radius*math.cos(math.pi/6)/500. * c_l[2]

        if x1 > 140 and x1 < 150:
            l_intercept.append(y1)
        
        ##### ##### ##### ##### ##### #####
        # Triangle area
        nvert = 3
        l_vertx = [[0,0,-150],[0,-150,-150],[0,-150,0],[0,0,150],[0,150,150],[0,150,0]]
        l_verty = [[0,-150,-75],[0,-75,75],[0,75,150],[0,150,75],[0,75,-75],[0,-75,-150]]
        tri = 0
        for vertx,verty in zip(l_vertx,l_verty):
            if func.isPointInPoly(nvert, vertx, verty, x1, y1):# or  func.isPointInPoly(nvert, vertx, verty, x2, y2):
                h_area_density_1.Fill(tri)
            tri += 1
            
        # Inner hexagon area
        nvert = 6
        l_vertx = [[0,-100,-100,0,100,100]]
        l_verty = [[-100,-50,50,100,50,-50]]
        for vertx,verty in zip(l_vertx,l_verty):
            if func.isPointInPoly(nvert, vertx, verty, x1, y1):# or  func.isPointInPoly(nvert, vertx, verty, x2, y2):
                h_area_density_2.Fill(1)
            else:
                h_area_density_2.Fill(0)

        # 2D track length
        h_2D_track_length.Fill(math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)))
        
        # 3D track length
        h_3D_track_length.Fill(math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)))

        ##### ##### ##### ##### ##### #####

        (mp_e, tm_e) = ([x1,y1],1)
        (mp_l, tm_l) = ([x2,y2],1)

        ints_e = func.getHexagon(mp_e)
        ints_l = func.getHexagon(mp_l)
        
        ints_all = func.combineAllPoints(ints_e, ints_l)
        ints_all = func.convex_hull(ints_all)
        ints_all = np.array(ints_all)
                        
        allmps_e.append(mp_e)
        allmps_l.append(mp_l)
        polyregion_all.append(ints_all)
        
        if tog.showprogress == 1:            
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

    intercept = np.average(l_intercept)
    print 'intercept',intercept
    func.plotAllRegionsAxis([allmps_e, allmps_l], polyregion_all, prefix,h_angle.GetMean(), intercept)
                
    if tog.locallaptop == 0:
        c1 = ROOT.TCanvas("c1","",900,600)
        h_angle.Draw()
        c1.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c1.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')
        c1.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.C')
        c2 = ROOT.TCanvas("c2","",900,600)
        h_area_density_1.SetMinimum(0.0)
        h_area_density_1.SetMaximum(500.0)
        h_area_density_1.DrawNormalized('ep')
        c2.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_triangles_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c2.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_triangles_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')
        c2.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_triangles_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.C')
        c3 = ROOT.TCanvas("c3","",900,600)
        h_area_density_2.SetMinimum(0.0)
        h_area_density_2.DrawNormalized('ep')
        c3.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_innerhex_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c3.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_innerhex_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')
        c3.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_innerhex_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.C')
        c4 = ROOT.TCanvas("c4","",900,600)
        h_2D_track_length.Draw()
        c4.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c4.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')
        c4.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.C')
        c4.SetLogy()
        c4.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'_log.pdf')
        c4.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'_log.png')
        c5 = ROOT.TCanvas("c5","",900,600)
        h_3D_track_length.Draw()
        c5.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.pdf')
        c5.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.png')
        c5.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'.C')
        c5.SetLogy()
        c5.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'_log.pdf')
        c5.Print('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+prefix[:-1]+'_'+str(tog.badarealimit)+'_log.png')

    if tog.createlog == 1: func.writeLog_4(timestamp[0], len(subevent), n_events, np.max(triareas_all), np.average(triareas_all), np.median(triareas_all), len(badevents))
    
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
ec_beam = []
lc_beam = []
subevent_beam = []

ec_cosmics = []
lc_cosmics = []
subevent_cosmics = []

n_events = -666

if tog.locallaptop == 0:
    # Read from ROOT file
    if len(sys.argv) < 2:
        raise ValueError('Missing input file.')
        
    f = ROOT.TFile(sys.argv[1])
    t = ROOT.TTree()
    f.GetObject("tracks",t)
    os.system('mkdir -p plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/events')
    
    ev = 1
    n_events = t.GetEntries()
    i = 0
    
    for e in t:
        if e.GetListOfBranches().FindObject('first_area'):
            first_areas = e.first_areas
            last_areas = e.first_areas
        else:
            first_areas = [0 for x in e.beam]
            last_areas = [0 for x in e.beam]
        if e.GetListOfBranches().FindObject('first_hit_z'):
            first_hit_z = e.first_hit_z
            last_hit_z = e.last_hit_z
        else:
            first_hit_z = [0 for x in e.beam]
            last_hit_z = [0 for x in e.beam]

        for hx1,hy1,hz1,hx2,hy2,hz2,a1,a2,beam in zip(e.first_hit_x,e.first_hit_y,first_hit_z,e.last_hit_x,e.last_hit_y,last_hit_z,first_areas,last_areas,e.beam):
            if a1 < tog.badarealimit and a2 < tog.badarealimit:
                if beam:
                    ec_beam.append([int(hx1),int(hy1),int(hz1)])
                    lc_beam.append([int(hx2),int(hy2),int(hz2)])
                    subevent_beam.append(1.1)
                else:
                    ec_cosmics.append([int(hx1),int(hy1),int(hz1)])
                    lc_cosmics.append([int(hx2),int(hy2),int(hz2)])
                    subevent_cosmics.append(1.1)
else:
    #tracks_2D_12100.root 
    ew_beam = loc.ec_beam
    lw_beam = loc.lc_beam
    subevent_beam = loc.subevent_beam
    ew_cosmics = loc.ec_cosmics
    lw_cosmics = loc.lc_cosmics
    subevent_cosmics = loc.subevent_cosmics

#ew_beam = [[ 8,110,311],[306,186,325],[278,-999,239]]
#lw_beam = [[21,126,318],[319,208,318],[305,-999,263]]
#ew_beam = [[278,-999,239],[270,110,310],[306,150,325],[8,145,-999]]
#lw_beam = [[305,-999,263],[280,126,325],[319,185,318],[21,165,-999]]
#subevent_beam = [1.0,1.1,1.2,1.3,1.4,1.5,.16,1.7,1.8]
#ew_cosmics = []
#lw_cosmics = []
#subevent_cosmics = []
    
###############################################################################
# run
###############################################################################

timestamp = func.getTime()

if tog.createlog == 1:
    func.writeLog_1(timestamp, ec_beam, lc_beam, subevent_beam, ec_cosmics, lc_cosmics, subevent_cosmics)

ans1 = wiresMain(n_events, ec_beam,lc_beam,subevent_beam,'plot_beam_wires_')
ans2 = wiresMain(n_events, ec_cosmics,lc_cosmics,subevent_cosmics,'plot_cosmics_wires_')
