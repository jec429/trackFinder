# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 18:38:49 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

"""

import os,sys,glob
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

import math
from array import array


###############################################################################
# initial start/setup
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

    nfiles = 0
    t = ROOT.TChain('tracks')
    for fn in sys.argv:
        if '.root' in fn:
            nfiles += 1    
            f = ROOT.TFile(fn)
            t.Add(fn)
            
            
            
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
    inputfile = "mumbojumbo/test.root"
    ec_beam = loc.ec_beam
    lc_beam = loc.lc_beam
    subevent_beam = loc.subevent_beam
    ec_cosmics = loc.ec_cosmics
    lc_cosmics = loc.lc_cosmics
    subevent_cosmics = loc.subevent_cosmics

if nfiles == 1:
    inputfile = sys.argv[1].split('/')[-1]
else:
    inputfile = 'tree_combined'

###############################################################################

os.system('mkdir -p plots/'+inputfile.replace('.root','')+'/limit'+str(tog.badarealimit))

timestamp = func.getTime()

func.writeWireID_list(timestamp, inputfile, ec_beam, lc_beam, subevent_beam, ec_cosmics, lc_cosmics, subevent_cosmics)

###############################################################################

def wiresRecon(n_events,ec,lc,subevent,ts,inf,pf):
    alltracks = allTrackInfo(ts, pf, inf, -1, n_events)
    func.writeAllTrackInfo_1(alltracks)

    x, y, z =  array( 'd' ), array( 'd' ), array( 'd' )

    if tog.locallaptop == 0:# and tog.plt_angles == 1:
        h_angle = ROOT.TH1F("h_angle","",100,-1.6,1.6)
        h_2D_track_length = ROOT.TH1F("h_2D_track_length","",100,0,1000)
        h_3D_track_length = ROOT.TH1F("h_3D_track_length","",100,0,1000)
        h_area_density_1 = ROOT.TH1F('h_area_density_1','',6,0,6)
        h_area_density_2 = ROOT.TH1F('h_area_density_2','',2,0,2)     

        h_dist_to_axis_x = ROOT.TH1F('h_dist_to_axis_x','',50,0,1000)     
        h_dist_to_axis_x2 = ROOT.TH1F('h_dist_to_axis_x2','',50,0,1000)     
        h_dist_to_axis_y = ROOT.TH1F('h_dist_to_axis_y','',100,-500,500)
        h_dist_to_axis_z = ROOT.TH1F('h_dist_to_axis_z','',100,-1000,300)
        h_dist_to_axis_z2 = ROOT.TH1F('h_dist_to_axis_z2','',100,-200,300)
        h_dist_to_axis = ROOT.TH2F('h_dist_to_axis','',50,0,1000,50,-500,500)
        h_dist_to_axis_xz = ROOT.TH2F('h_dist_to_axis_xz','',50,0,1000,50,-200,300)
        h_dist_to_axis_yz = ROOT.TH2F('h_dist_to_axis_yz','',50,-500,500,50,-200,300)
        h_dist_to_axis_xyz = ROOT.TH3F('h_dist_to_axis_xyz','',50,0,1000,50,-500,500,50,-200,300)
        
    
    for c_e,c_l,se in zip(ec,lc,subevent):
        newtrack = trackInfo(se=se)
        
        x1 = const.radius*const.dd/250. * c_e[0]
        y1 = const.radius*const.dd/250. * c_e[1]
        z1 = const.radius*const.dd/250. * c_e[2]
        x2 = const.radius*const.dd/250. * c_l[0]
        y2 = const.radius*const.dd/250. * c_l[1]
        z2 = const.radius*const.dd/250. * c_l[2]

        if 150 - x1 < 150 - x2:
            if y1 < 50 and y1 > -100 and x1 > 0:
                x.append(x1)
                y.append(y1)
                z.append(z1)
        else:
            if y2 < 50 and y2 > -100 and x2 > 0:
                x.append(x2)
                y.append(y2)
                z.append(z2)
        
        if x1 > 140 and x1 < 150:
            newtrack.li = y1
        
        ##### ##### ##### ##### ##### #####
        # Triangle area
        nvert = 3
        l_vertx = [[0,0,-150],[0,-150,-150],[0,-150,0],[0,0,150],[0,150,150],[0,150,0]]
        l_verty = [[0,-150,-75],[0,-75,75],[0,75,150],[0,150,75],[0,75,-75],[0,-75,-150]]
        tri = 0
        h_ad1s = []
        for vertx,verty in zip(l_vertx,l_verty):
            if func.isPointInPoly(nvert, vertx, verty, x1, y1):# or  func.isPointInPoly(nvert, vertx, verty, x2, y2):
                if tog.locallaptop == 0:
                    h_area_density_1.Fill(tri)
                else:
                    h_ad1s.append(tri)
            tri += 1
        newtrack.h_ad1 = h_ad1s
            
        # Inner hexagon area
        nvert = 6
        l_vertx = [[0,-100,-100,0,100,100]]
        l_verty = [[-100,-50,50,100,50,-50]]
        h_ad2s = []
        for vertx,verty in zip(l_vertx,l_verty):
            if func.isPointInPoly(nvert, vertx, verty, x1, y1):# or  func.isPointInPoly(nvert, vertx, verty, x2, y2):
                if tog.locallaptop == 0:
                    h_area_density_2.Fill(1)
                else:
                    h_ad2s.append(1)
            else:
                if tog.locallaptop == 0:
                    h_area_density_2.Fill(0)
                else:
                    h_ad2s.append(0)
        newtrack.h_ad2 = h_ad1s


        # 2D track length
        h2Dtl = math.sqrt((c_l[0]-c_e[0])*(c_l[0]-c_e[0])+(c_l[1]-c_e[1])*(c_l[1]-c_e[1]))
        # 3D track length
        h3Dtl = math.sqrt((c_l[0]-c_e[0])*(c_l[0]-c_e[0])+(c_l[1]-c_e[1])*(c_l[1]-c_e[1])+(c_l[2]-c_e[2])*(c_l[2]-c_e[2]))
        
        if tog.locallaptop == 0:
            # 2D track length
            h_2D_track_length.Fill(h2Dtl)
            # 3D track length
            h_3D_track_length.Fill(h3Dtl)
        else:
            newtrack.h2Dtl = h2Dtl
            newtrack.h3Dtl = h3Dtl

        ##### ##### ##### ##### ##### #####
        
        #WAIT! QUESTION! WHAT IS tm_e HERE?
        (mp_e, tm_e) = ([x1,y1],1)
        (mp_l, tm_l) = ([x2,y2],1)

        newtrack.mp_e = mp_e
        newtrack.mp_l = mp_l
        #newtrack.mt_e = sum(tm_e)
        #newtrack.mt_l = sum(tm_l)
        newtrack.mt_e = tm_e
        newtrack.mt_l = tm_l
        
        ang = func.getTrackAngle(mp_e,mp_l)
        newtrack.ang = ang
        if ang != -666 and tog.locallaptop == 0 and tog.plt_angles == 1:
            h_angle.Fill(ang)
        
        ints_e = func.getHexagon(mp_e)
        ints_l = func.getHexagon(mp_l)
        
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
                
        #if newtrack.ia_e < 0.433 or newtrack.ia_l < 0.433:
        alltracks.add_track(newtrack)
        
    if tog.locallaptop == 0 and tog.plt_angles == 1:
        ROOT.gStyle.SetOptStat(111111)
        c1 = ROOT.TCanvas("c1","",900,600)
        h_angle.Draw()
        h_angle.GetXaxis().SetTitle("track angle w.r.t x-axis (rad)")
        plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_angle_'+pf+'_'+str(tog.badarealimit)+'_reco'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c1.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c1.Print(plotname+'.png')
        c1.Print(plotname+'.C')

        c2 = ROOT.TCanvas("c2","",900,600)
        h_area_density_1.SetMinimum(0.0)
        #h_area_density_1.SetMaximum(500.0)
        h_area_density_1.DrawNormalized('ep')
        h_area_density_1.GetXaxis().SetTitle('Triangle index')
        h_area_density_1.GetYaxis().SetTitle('Relative density')
        plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_triangles_'+pf+'_'+str(tog.badarealimit)+'_reco'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c2.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c2.Print(plotname+'.png')
        c2.Print(plotname+'.C')

        c3 = ROOT.TCanvas("c3","",900,600)
        h_area_density_2.SetMinimum(0.0)
        h_area_density_2.DrawNormalized('ep')
        h_area_density_2.GetXaxis().SetTitle('Inside inner hexagon')
        h_area_density_2.GetYaxis().SetTitle('Relative density')
        plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_area_density_innerhex_'+pf+'_'+str(tog.badarealimit)+'_reco'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c3.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c3.Print(plotname+'.png')
        c3.Print(plotname+'.C')

        c4 = ROOT.TCanvas("c4","",900,600)
        h_2D_track_length.Draw()
        h_2D_track_length.GetXaxis().SetTitle('track length 2D [mm]')
        plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_2D_track_length_'+pf+'_'+str(tog.badarealimit)+'_reco'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c4.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c4.Print(plotname+'.png')
        c4.Print(plotname+'.C')
        c4.SetLogy()
        plotname = plotname + '_log'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c4.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c4.Print(plotname+'.png')

        c5 = ROOT.TCanvas("c5","",900,600)
        h_3D_track_length.Draw()
        h_3D_track_length.GetXaxis().SetTitle('track length 3D [mm]')
        plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_3D_track_length_'+pf+'_'+str(tog.badarealimit)+'_reco'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c5.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c5.Print(plotname+'.png')
        c5.Print(plotname+'.C')
        c5.SetLogy()
        plotname = plotname + '_log'
        if tog.saveplotformat == 1 or tog.saveplotformat == 2: c5.Print(plotname+'.pdf')
        if tog.saveplotformat == 0 or tog.saveplotformat == 2: c5.Print(plotname+'.png')
    
    # QUESTION! HOW INTEGRATE
    #func.writeAllTrackInfo_3(alltracks, -1)
    
    #func.writeTrackInfo(alltracks)
    gr = ROOT.TGraph(len(x),x,y)
    gr.Fit('pol1','FQ')
    f1 = gr.GetFunction('pol1')
    new_coords = func.calculateNewCoordinates(alltracks,f1.GetParameter(1),f1.GetParameter(0))
    for ep,lp,c,beamz in zip(alltracks.mp_e,alltracks.mp_l,new_coords,z):

        beamx = c[0] if 150 - ep[0] < 150 - lp[0] else c[2]
        beamy = c[1] if 150 - ep[0] < 150 - lp[0] else c[3]
        if abs(beamz + 235) == min([abs(beamz + 235), abs(beamz + 142), abs(beamz + 49)]):
            beamz_c = beamz + 235
        elif abs(beamz + 142) == min([abs(beamz + 235), abs(beamz + 142), abs(beamz + 49)]):
            beamz_c = beamz + 142
        elif abs(beamz + 49) == min([abs(beamz + 235), abs(beamz + 142), abs(beamz + 49)]):
            beamz_c = beamz + 49
        else:
            beamz_c = -999

        #beamx = c[2]
        #beamy = c[3]
        
        h_dist_to_axis_x.Fill((150-beamx)/(const.radius*const.dd/250.))
        if abs(beamy) < 50 and abs(beamz_c) < 50:
            h_dist_to_axis_x2.Fill((150-beamx)/(const.radius*const.dd/250.))
        h_dist_to_axis_y.Fill(beamy/(const.radius*const.dd/250.))
        h_dist_to_axis_z.Fill(beamz/(const.radius*const.dd/250.))
        h_dist_to_axis_z2.Fill(beamz_c/(const.radius*const.dd/250.))
        h_dist_to_axis.Fill((150-beamx)/(const.radius*const.dd/250.),beamy/(const.radius*const.dd/250.))
        h_dist_to_axis_xz.Fill((150-beamx)/(const.radius*const.dd/250.),beamz_c/(const.radius*const.dd/250.))
        h_dist_to_axis_yz.Fill(beamy/(const.radius*const.dd/250.),beamz_c/(const.radius*const.dd/250.))
        h_dist_to_axis_xyz.Fill((150-beamx)/(const.radius*const.dd/250.),beamy/(const.radius*const.dd/250.),beamz_c/(const.radius*const.dd/250.))
        
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_x.Draw()
    h_dist_to_axis_x.GetXaxis().SetTitle('x distance to axis [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_x_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')
    c6.SetLogy()
    plotname = plotname + '_log'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    
    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_z.Draw()
    g1 = ROOT.TF1("g1","[2]/(3.14*[0])*([0]*[0]/((x-[1])*(x-[1])+[0]*[0]))",-850,-750) # Lorentz/Cauchy distribution
    g1.SetParameter(0,10.0)
    g1.SetParameter(1,-800.0)
    g1.SetParameter(2,10000.0)
    h_dist_to_axis_z.Fit('g1','R')
    g2 = ROOT.TF1("g2","[2]/(3.14*[0])*([0]*[0]/((x-[1])*(x-[1])+[0]*[0]))",-530,-430) # Lorentz/Cauchy distribution
    g2.SetParameter(0,10.0)
    g2.SetParameter(1,-500.0)
    g2.SetParameter(2,10000.0)
    h_dist_to_axis_z.Fit('g2','R')
    g3 = ROOT.TF1("g3","[2]/(3.14*[0])*([0]*[0]/((x-[1])*(x-[1])+[0]*[0]))",-210,-110) # Lorentz/Cauchy distribution
    g3.SetParameter(0,10.0)
    g3.SetParameter(1,-160.0)
    g3.SetParameter(2,10000.0)
    h_dist_to_axis_z.Fit('g3','R')
    g1.Draw('same')
    g2.Draw('same')
    g3.Draw('same')
    h_dist_to_axis_z.GetXaxis().SetTitle('z position [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_z_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')
    c6.SetLogy()
    plotname = plotname + '_log'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_z2.Draw()
    g1 = ROOT.TF1("g1","[2]/(3.14*[0])*([0]*[0]/((x-[1])*(x-[1])+[0]*[0]))",-40,40) # Lorentz/Cauchy distribution
    g1.SetParameter(0,10.0)
    g1.SetParameter(1,1.0)
    g1.SetParameter(2,10000.0)
    h_dist_to_axis_z2.Fit('g1','R')
    h_dist_to_axis_z2.GetXaxis().SetTitle('z position [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_z2_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')
    c6.SetLogy()
    plotname = plotname + '_log'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_x2.Draw("ep")
    h_dist_to_axis_x2.Fit("expo")
    h_dist_to_axis_x2.GetXaxis().SetTitle('x distance to axis [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_x2_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')
    c6.SetLogy()
    plotname = plotname + '_log'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_y.Draw()   
    #f1 = ROOT.TF1("f1","[0]*14.0/22.0*[1]*[1]*[2]*[2]/(((x*x) - ([2]*[2]))*((x*x) - ([2]*[2])) + x*x*x*x*(([1]*[1])/([2]*[2])))",-50,50) # Breit-Wigner
    f1 = ROOT.TF1("f1","[2]/(3.14*[0])*([0]*[0]/((x-[1])*(x-[1])+[0]*[0]))",-50,50) # Lorentz/Cauchy distribution
    #f1 = ROOT.TF1("f1","gaus",-50,50)
    f1.SetParameter(0,50.0)
    f1.SetParameter(1,1.0)
    f1.SetParameter(2,10000.0)
    h_dist_to_axis_y.Fit('f1','R')
    h_dist_to_axis_y.GetXaxis().SetTitle('y distance to axis [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_y_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')
    c6.SetLogy()
    plotname = plotname + '_log'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    
    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,900)
    ROOT.gStyle.SetOptStat(0)
    h_dist_to_axis.Draw('colz')
    h_dist_to_axis.GetXaxis().SetTitle('x distance to axis [mm]')
    h_dist_to_axis.GetYaxis().SetTitle('y distance to axis [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_xz.Draw('colz')
    h_dist_to_axis_xz.GetXaxis().SetTitle('x position [mm]')
    h_dist_to_axis_xz.GetYaxis().SetTitle('z position [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_xz_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_yz.Draw('colz')
    h_dist_to_axis_yz.GetXaxis().SetTitle('y position [mm]')
    h_dist_to_axis_yz.GetYaxis().SetTitle('z position [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_yz_'+pf+'_'+str(tog.badarealimit)+'_reco'
    if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')

    c6 = 0
    c6 = ROOT.TCanvas("c6","",900,600)
    h_dist_to_axis_xyz.Draw("BOX2Z")
    #h_dist_to_axis_yz.GetXaxis().SetTitle('y position [mm]')
    #h_dist_to_axis_yz.GetYaxis().SetTitle('z position [mm]')
    plotname = 'plots/'+inf.replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+'h_dist_to_axis_xyz_'+pf+'_'+str(tog.badarealimit)+'_reco'
    #if tog.saveplotformat == 1 or tog.saveplotformat == 2: c6.Print(plotname+'.pdf')
    #if tog.saveplotformat == 0 or tog.saveplotformat == 2: c6.Print(plotname+'.png')
    c6.Print(plotname+'.C')

    
    #if tog.plt_lines == 1 or tog.plt_shading == 1:
        #func.plotAllRegions(alltracks)
        #func.plotAllRegionsAxis_newFormat(alltracks,f1.GetParameter(1),f1.GetParameter(0))
    #return alltracks

#

###############################################################################
# run
###############################################################################

tracks_beam = wiresRecon(n_events, ec_beam,lc_beam,subevent_beam,timestamp,inputfile,'beam')
#tracks_cosmics = wiresRecon(n_events, ec_cosmics,lc_cosmics,subevent_cosmics,timestamp,inputfile,'cosmics')
