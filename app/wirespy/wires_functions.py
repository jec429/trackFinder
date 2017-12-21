# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 13:20:44 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: functions that wires.py uses
"""

from datetime import datetime
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#from pathlib import Path

import wires_constants as const
import wires_toggle as tog


###############################################################################
# classes
###############################################################################

class trackInfo:
    def __init__(self, subevent, wireID_e, wireID_l, midpointORintersection_e, midpointORintersection_l, wirematchtype_e, wirematchtype_l, intersectioncoords_e, intersectioncoords_l, intersectionarea_e, intersectionarea_l, polygonarea, longesttracklength, trackangle, h_2D_track_length, h_3D_track_length, h_area_density_1, h_area_density_2):
        self.se = subevent
        self.w_e = wireID_e
        self.w_l = wireID_l
        self.mp_e = midpointORintersection_e
        self.mp_l = midpointORintersection_l
        self.mt_e = wirematchtype_e
        self.mt_l = wirematchtype_l
        self.ic_e = intersectioncoords_e
        self.ic_l = intersectioncoords_l
        self.ia_e = intersectionarea_e
        self.ia_l = intersectionarea_l
        self.pa = polygonarea
        self.tl = longesttracklength
        self.ang = trackangle
        self.li = l_intercept
        self.h2Dtl = h_2D_track_length
        self.h3Dtl = h_3D_track_length
        self.h_ad1 = h_area_density_1
        self.h_ad2 = h_area_density_2

class allTrackInfo:
    def __init__(self, timestamp, prefix, inputfile, limit, num_events):
        self.ts = timestamp
        self.pf = prefix
        self.inf = inputfile
        self.bal = limit
        self.ne = num_events
        self.se = []
        self.w_e = []
        self.w_l = []
        self.mp_e = []
        self.mp_l = []
        self.mt_e = []
        self.mt_l = []
        self.ic_e = []
        self.ic_l = []
        self.ia_e = []
        self.ia_l = []
        self.pa = []
        self.tl = []
        self.ang = []
        self.li = []
        self.h2Dtl = []
        self.h3Dtl = []
        self.h_ad1 = []
        self.h_ad2 = []
        
    def add_track(self, trackInfo):
        self.se.append(trackInfo.se)
        self.w_e.append(trackInfo.w_e)
        self.w_l.append(trackInfo.w_l)
        self.mp_e.append(trackInfo.mp_e)
        self.mp_l.append(trackInfo.mp_l)
        self.mt_e.append(trackInfo.mt_e)
        self.mt_l.append(trackInfo.mt_l)
        self.ic_e.append(trackInfo.ic_e)
        self.ic_l.append(trackInfo.ic_l)
        self.ia_e.append(trackInfo.ia_e)
        self.ia_l.append(trackInfo.ia_l)
        self.pa.append(trackInfo.pa)
        self.tl.append(trackInfo.tl)
        self.ang.append(trackInfo.ang)
        self.li.append(trackInfo.li)
        self.h2Dtl.append(trackInfo.h2Dtl)
        self.h3Dtl.append(trackInfo.h3Dtl)
        self.h_ad1.append(trackInfo.h_ad1)
        self.h_ad2.append(trackInfo.h_ad2)

    def print_basics(self):
        print("timestamp =\t" + str(self.ts))
        print("prefix =\t" + str(self.pf))
        print("inputfile =\t" + str(self.inf))
        print("limit =\t\t" + str(self.bal))
        print("num_events =\t" + str(self.ne))
        print("num of tracks =\t" + str(len(self.se)))

###############################################################################
# functions
###############################################################################

############################################################ plotting functions

def plotAxis():
    color = (200, 200, 200)
    color = '#%02x%02x%02x' % color
    plt.plot(const.xaxis[0], const.xaxis[1], color, ls='-')
    plt.plot(const.uaxis[0], const.uaxis[1], color, ls='-')
    plt.plot(const.vaxis[0], const.vaxis[1], color, ls='-')

#

def plotMPs(mp_e, mp_l):
    plt.plot(mp_e[0], mp_e[1], const.startcolor, marker='.', markersize=const.ms)
    plt.plot(mp_l[0], mp_l[1], const.endcolor, marker='.', markersize=const.ms)

#

def plotAllRegions(alltracks):
    #currently have fixed stylistic settings
    mps_e = np.array(alltracks.mp_e)
    mps_l = np.array(alltracks.mp_l)
    polyregions = alltracks.pa
    
    if tog.plt_shading == 1:
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
    
    if tog.plt_lines == 1:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, aspect='equal')
        ax2.add_patch(
            patches.RegularPolygon(
                (0, 0),      # (x,y)
                6,           # number of vertices
                const.radius,      # radius
                fill = False
            )
        )
        ax2.set_xlim(-200,200)
        ax2.set_ylim(-200,200)
    
    if const.plotaxis == 1: plotAxis()
    
    if tog.plt_shading == 1:
        for aa in range(len(polyregions)):
            ints_all = polyregions[aa]
            ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=const.hotshading, lw=0))

    #if tog.plt_lines == 1:
        #plt.plot([mps_e[:,0], mps_l[:,0]], [mps_e[:,1], mps_l[:,1]], const.midcolor, ls='-', lw=const.lw)
        #plt.plot(mps_e[:,0], mps_e[:,1], const.startcolor, marker='.', markersize=const.ms, lw=0)
        #plt.plot(mps_l[:,0], mps_l[:,1], const.endcolor, marker='.', markersize=const.ms, lw=0)
    
    inf = alltracks.inf
    filename = 'plots/'+inf.replace('.root','').split('/')[-1]+'/'
    if alltracks.bal != -1:
        filename = filename + "/limit" + str(alltracks.bal) + "/limit" + str(alltracks.bal) + "_"
    filename = filename + alltracks.pf + '_'
    
    if tog.plt_lines == 1:
        if tog.saveplotformat == 0 or tog.saveplotformat == 2:
            plt.savefig(filename + 'lines_midpoints' + '.png', bbox_inches='tight')
        if tog.saveplotformat == 1 or tog.saveplotformat == 2:
            plt.savefig(filename + 'lines_midpoints' + '.pdf', bbox_inches='tight')
        plt.close()
    
    if tog.plt_shading == 1:
        #plt.plot(mps_e[:,0], mps_e[:,1], const.startcolor, marker='.', markersize=const.ms, lw=0)
        #plt.plot(mps_l[:,0], mps_l[:,1], const.endcolor, marker='.', markersize=const.ms, lw=0)
        
        if tog.saveplotformat == 0 or tog.saveplotformat == 2:
            plt.savefig(filename + 'shading_midpoints'+'.png', bbox_inches='tight')
        if tog.saveplotformat == 1 or tog.saveplotformat == 2:
            plt.savefig(filename + 'shading_midpoints'+'.pdf', bbox_inches='tight')
        plt.close()

#

def plotAllRegionsAxis(allmps, polyregion_all, prefix, beam_angle, intercept):
    import sys
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
    if const.plotaxis == 1: plotAxis()
    allmps_e = allmps[0]
    allmps_l = allmps[1]
    if tog.plt_shading == 1:
        for aa in range(len(allmps_e)):
            ints_all = polyregion_all[aa]
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=const.hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
    if tog.plt_lines == 1:
        for aa in range(len(allmps_e)):
            mp_e = allmps_e[aa]
            mp_l = allmps_l[aa]
            plt.plot([mp_e[0], mp_l[0]], [mp_e[1], mp_l[1]], const.midcolor, ls='-', lw=const.lw)
    if tog.plt_dots == 1:
        for aa in range(len(allmps_e)):
            mp_e = allmps_e[aa]
            mp_l = allmps_l[aa]
            plotMPs(mp_e, mp_l)

    b = intercept-math.tan(beam_angle)*(150)
    y1 = math.tan(beam_angle)*(-150)+b
    y2 = math.tan(beam_angle)*(150)+b
    if tog.locallaptop == 0:
        print('axis',math.tan(beam_angle),b,y1,y2)
    plt.plot([-150,150],[y1,y2],linewidth=2.0)

    if tog.locallaptop == 1:
        plt.savefig(prefix+str(tog.badarealimit)+'_'+str(tog.plt_dots)+str(tog.plt_lines)+str(tog.plt_shading)+'_midpoints_axis.pdf', bbox_inches='tight')
    else:
        plt.savefig('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+prefix+str(tog.badarealimit)+'_'+str(tog.plt_dots)+str(tog.plt_lines)+str(tog.plt_shading)+'_midpoints_axis.pdf', bbox_inches='tight')
        plt.savefig('plots/'+sys.argv[1].replace('.root','').split('/')[-1]+'/limit'+str(tog.badarealimit)+'/'+prefix+str(tog.badarealimit)+'_'+str(tog.plt_dots)+str(tog.plt_lines)+str(tog.plt_shading)+'_midpoints_axis.png', bbox_inches='tight')

#

def plotAllRegionsAxis_newFormat(alltracks, beam_angle, intercept):
    mps_e = np.array(alltracks.mp_e)
    mps_l = np.array(alltracks.mp_l)
    polyregions = alltracks.pa
    
    b = intercept-math.tan(beam_angle)*(150)
    y1 = math.tan(beam_angle)*(-150)+b
    y2 = math.tan(beam_angle)*(150)+b
    if tog.locallaptop == 0:
        print('axis',math.tan(beam_angle),b,y1,y2)

    if tog.plt_shading == 1:
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
    
    if tog.plt_lines == 1:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, aspect='equal')
        ax2.add_patch(
            patches.RegularPolygon(
                (0, 0),      # (x,y)
                6,           # number of vertices
                const.radius,      # radius
                fill = False
            )
        )
        ax2.set_xlim(-200,200)
        ax2.set_ylim(-200,200)
    
    if const.plotaxis == 1: plotAxis()
    
    if tog.plt_shading == 1:
        for aa in range(len(polyregions)):
            ints_all = polyregions[aa]
            ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=const.hotshading, lw=0))

    if tog.plt_lines == 1:
        plt.plot([mps_e[:,0], mps_l[:,0]], [mps_e[:,1], mps_l[:,1]], const.midcolor, ls='-', lw=const.lw)
        plt.plot(mps_e[:,0], mps_e[:,1], const.startcolor, marker='.', markersize=const.ms, lw=0)
        plt.plot(mps_l[:,0], mps_l[:,1], const.endcolor, marker='.', markersize=const.ms, lw=0)
        plt.plot([-150,150],[y1,y2],linewidth=2.0)
    
    inf = alltracks.inf
    filename = 'plots/'+inf.replace('.root','').split('/')[-1]+'/'
    if alltracks.bal != -1:
        filename = filename + "/limit" + str(alltracks.bal) + "/limit" + str(alltracks.bal) + "_"
    filename = filename + alltracks.pf + '_'
    
    if tog.plt_lines == 1:
        if tog.saveplotformat == 0 or tog.saveplotformat == 2:
            plt.savefig(filename + 'lines_midpoints_axis' + '.png', bbox_inches='tight')
        if tog.saveplotformat == 1 or tog.saveplotformat == 2:
            plt.savefig(filename + 'lines_midpoints_axis' + '.pdf', bbox_inches='tight')
        plt.close()
    
    if tog.plt_shading == 1:
        plt.plot(mps_e[:,0], mps_e[:,1], const.startcolor, marker='.', markersize=const.ms, lw=0)
        plt.plot(mps_l[:,0], mps_l[:,1], const.endcolor, marker='.', markersize=const.ms, lw=0)
        plt.plot([-150,150],[y1,y2],linewidth=2.0)

        if tog.saveplotformat == 0 or tog.saveplotformat == 2:
            plt.savefig(filename + 'shading_midpoints_axis'+'.png', bbox_inches='tight')
        if tog.saveplotformat == 1 or tog.saveplotformat == 2:
            plt.savefig(filename + 'shading_midpoints_axis'+'.pdf', bbox_inches='tight')
        plt.close()

#################################################################### create log

#def getTime(inputfile):
def getTime():
    s_y = datetime.now().strftime('%y')
    s_m = datetime.now().strftime('%m')
    s_d = datetime.now().strftime('%d')
    s_H = datetime.now().strftime('%H')
    s_M = datetime.now().strftime('%M')
    s_S = datetime.now().strftime('%S')
    filename = s_y + s_m + s_d + s_H + s_M + s_S
    ts = s_d + "/" + s_m + "/" + s_y + " " + s_H + ":" + s_M + ":" + s_S
    return [filename, ts]

#

def writeWireID_list(ts, inf, ew_beam, lw_beam, subevent_beam, ew_cosmics, lw_cosmics, subevent_cosmics):
    filename = 'plots/' + inf.replace('.root','').split('/')[-1] + "/" + inf.replace('.root','_').split('/')[-1] + "wireIDlist.txt"
    #print(filename)
    text_file = open(filename, "w")
    text_file.write("ew_beam = " + str(ew_beam) + "\n")
    text_file.write("lw_beam = " + str(lw_beam) + "\n")
    text_file.write("subevent_beam = " + str(subevent_beam) + "\n\n\n\n")
    text_file.write("ew_cosmics = " + str(ew_cosmics) + "\n")
    text_file.write("lw_cosmics = " + str(lw_cosmics) + "\n")
    text_file.write("subevent_cosmics = " + str(subevent_cosmics) + "\n")
    text_file.close()

#

def writeAllTrackInfo_1(allTrackInfo):
    inf = allTrackInfo.inf
    filename = 'plots/' + inf.replace('.root','').split('/')[-1] 
    if allTrackInfo.bal != -1:
        filename = filename + "/limit" + str(allTrackInfo.bal) + "/"
    filename = filename + "/log_" + allTrackInfo.ts[0] + "_"  + inf.replace('.root','_').split('/')[-1] + allTrackInfo.pf
    if allTrackInfo.bal == -1:
        filename = filename + "_allTrackInfo.txt"
    else:
        filename = filename + "_limit" + str(allTrackInfo.bal) + "_trackInfo.txt"
    text_file = open(filename, "w")
    text_file.write("timestamp: " + allTrackInfo.ts[1] + "\n")
    text_file.write("file: " + inf + "\n")
    text_file.write(allTrackInfo.pf + "\n")
    #text_file.write("subevent\tearly wire ID [x,u,v]\tlate wire ID [x,u,v]\tearly midpoint/intersection\tlate midpoint/intersection\ttype of wire match (early)\ttype of wire match (late)\tearly intersection area\tlate intersection area\tpolygon area\tlongest track length\ttrack angle\t\n")
    text_file.close()

#
"""
def writeAllTrackInfo_2(track, inf, ts, pf):
    filename = 'plots/' + inf.replace('.root','').split('/')[-1] + "/log_" + ts[0] + "_" + inf.replace('.root','_').split('/')[-1] + pf + "_allTrackInfo.txt"
    text_file = open(filename, "a")
    text_file.write("\n" + str(track.se) + "\t" + str(track.w_e) + "\t" + str(track.w_l) + "\t")
    text_file.write(str(track.mp_e) + "\t" + str(track.mp_l) + "\t" + str(track.mt_e) + "\t" + str(track.mt_l) + "\t")
    text_file.write(str(track.ia_e) + "\t" + str(track.ia_l) + "\t" + str(track.pa) + "\t")
    text_file.write(str(track.tl) + "\t" + str(track.ang) + "\t")
    text_file.close()
"""
#

def writeAllTrackInfo_3(allTrackInfo, n_alltracks):
    inf = allTrackInfo.inf
    filename = 'plots/' + inf.replace('.root','').split('/')[-1] 
    if allTrackInfo.bal != -1:
        filename = filename + "/limit" + str(allTrackInfo.bal) + "/"
    filename = filename + "/log_" + allTrackInfo.ts[0] + "_"  + inf.replace('.root','_').split('/')[-1] + allTrackInfo.pf
    if allTrackInfo.bal == -1:
        filename = filename + "_allTrackInfo.txt"
    else:
        filename = filename + "_limit" + str(allTrackInfo.bal) + "_trackInfo.txt"
    text_file = open(filename, "a")
    text_file.write("\nnumber of tracks = " + str(len(allTrackInfo.se)))
    text_file.write("\nnumber of events = " + str(allTrackInfo.ne))
    text_file.write("\nnumber of tracks per event = " + str(len(allTrackInfo.se)/float(allTrackInfo.ne)))
    areas = allTrackInfo.ia_e + allTrackInfo.ia_l
    text_file.write("\n\nmax area = " + str(np.max(areas)))
    text_file.write("\nmin area = " + str(np.min(areas)))
    text_file.write("\naverage area = " + str(np.average(areas)))
    text_file.write("\nmedian area = " + str(np.median(areas)))
    text_file.write("\n\nmax track length = " + str(np.max(allTrackInfo.tl)))
    text_file.write("\naverage track length = " + str(np.average(allTrackInfo.tl)))
    text_file.write("\nmedian track length = " + str(np.median(allTrackInfo.tl)))
    if allTrackInfo.bal != -1:
        text_file.write("\n\n# of bad tracks (area(s) > " + str(allTrackInfo.bal) + ") = " + str(n_alltracks) + " - " + str(len(allTrackInfo.se)) + " = " + str(n_alltracks - len(allTrackInfo.se)))
    text_file.write("\n")
    text_file.close()

################################################# saving & opening allTrackInfo

def writeTrackInfo(allTrackInfo):
    inf = allTrackInfo.inf
    filename = 'plots/' + inf.replace('.root','').split('/')[-1]
    if allTrackInfo.bal == -1:
        filename = filename + "/allTrackInfo_"
    else:
        filename = filename + "/limit" + str(allTrackInfo.bal) + "/limit" + str(allTrackInfo.bal) + "_trackInfo_"
    filename = filename + inf.replace('.root','_').split('/')[-1] + allTrackInfo.pf + ".txt"
    text_file = open(filename, "w")
    text_file.write(str(allTrackInfo.ts) + "\n")
    text_file.write(str(allTrackInfo.pf) + "\n")
    text_file.write(str(allTrackInfo.inf) + "\n")
    text_file.write(str(allTrackInfo.bal) + "\n")
    text_file.write(str(allTrackInfo.ne) + "\n")
    text_file.write(str(allTrackInfo.se) + "\n")
    text_file.write(str(allTrackInfo.w_e) + "\n")
    text_file.write(str(allTrackInfo.w_l) + "\n")
    text_file.write(str(allTrackInfo.mp_e) + "\n")
    text_file.write(str(allTrackInfo.mp_l) + "\n")
    text_file.write(str(allTrackInfo.mt_e) + "\n")
    text_file.write(str(allTrackInfo.mt_l) + "\n")
    text_file.write(str(allTrackInfo.ic_e) + "\n")
    text_file.write(str(allTrackInfo.ic_l) + "\n")
    text_file.write(str(allTrackInfo.ia_e) + "\n")
    text_file.write(str(allTrackInfo.ia_l) + "\n")
    text_file.write(str(allTrackInfo.pa) + "\n")
    text_file.write(str(allTrackInfo.tl) + "\n")
    text_file.write(str(allTrackInfo.ang) + "\n")
    text_file.close()

def readAllTrackInfo(inf, pf):
    filename = 'plots/' + inf.replace('.root','').split('/')[-1] + "/allTrackInfo_" + inf.replace('.root','_').split('/')[-1] + pf + ".txt"
    #file = Path(filename)
    #if file.exists():
    if 1 == 1:
        text_file = open(filename,"r")
        timestamp = eval(text_file.readline())
        prefix = text_file.readline()
        prefix = prefix.partition("\n")[0]
        inputfile = text_file.readline()
        inputfile = inputfile.partition("\n")[0]
        limit = eval(text_file.readline())
        num_events = eval(text_file.readline())
        allTracks = allTrackInfo(timestamp, prefix, inputfile, limit, num_events)
        allTracks.se = eval(text_file.readline())
        allTracks.w_e = eval(text_file.readline())
        allTracks.w_l = eval(text_file.readline())
        allTracks.mp_e = eval(text_file.readline())
        allTracks.mp_l = eval(text_file.readline())
        allTracks.mt_e = eval(text_file.readline())
        allTracks.mt_l = eval(text_file.readline())
        allTracks.ic_e = eval(text_file.readline())
        allTracks.ic_l = eval(text_file.readline())
        allTracks.ia_e = eval(text_file.readline())
        allTracks.ia_l = eval(text_file.readline())
        allTracks.pa = eval(text_file.readline())
        allTracks.tl = eval(text_file.readline())
        allTracks.ang = eval(text_file.readline())
        text_file.close()
        return allTracks
    else:
        print("UH OH!!!! BAD THINGS! NO FILE EXISTS!")


############################################################ find intersections
"""
def getIntersectionUV(uu, vv):
    ux = uu[0]
    uy = uu[1]
    vx = vv[0]
    vy = vv[1]
    #ux1 = ux[0]
    #ux2 = ux[1]
    #uy1 = uy[0]
    #uy2 = uy[1]
    m_u = (uy[1] - uy[0]) / (ux[1] - ux[0])
    #vx1 = vx[0]
    #vx2 = vx[1]
    #vy1 = vy[0]
    #vy2 = vy[1]
    m_v = (vy[1] - vy[0]) / (vx[1] - vx[0])
    i_x = (uy[1] - m_u * ux[1] - vy[1] + m_v * vx[1]) / (m_v - m_u)
    i_y = vy[1] + m_v * (i_x - vx[1])
    return [i_x, i_y]
"""
#

def getIntersectionUV(uu, vv):
    ux = uu[0]
    uy = uu[1]
    vx = vv[0]
    vy = vv[1]
    #m_u = ((uy[1] - uy[0]) / (ux[1] - ux[0]))
    #m_v = ((vy[1] - vy[0]) / (vx[1] - vx[0]))
    i_x = (uy[1] - ((uy[1] - uy[0]) / (ux[1] - ux[0])) * ux[1] - vy[1] + ((vy[1] - vy[0]) / (vx[1] - vx[0])) * vx[1]) / (((vy[1] - vy[0]) / (vx[1] - vx[0])) - ((uy[1] - uy[0]) / (ux[1] - ux[0])))
    i_y = vy[1] + ((vy[1] - vy[0]) / (vx[1] - vx[0])) * (((uy[1] - ((uy[1] - uy[0]) / (ux[1] - ux[0])) * ux[1] - vy[1] + ((vy[1] - vy[0]) / (vx[1] - vx[0])) * vx[1]) / (((vy[1] - vy[0]) / (vx[1] - vx[0])) - ((uy[1] - uy[0]) / (ux[1] - ux[0])))) - vx[1])
    #i_y = vy[1] + ((vy[1] - vy[0]) / (vx[1] - vx[0])) * (i_x - vx[1])
    #i_x = (uy[1] - m_u * ux[1] - vy[1] + m_v * vx[1]) / (m_v - m_u)
    #i_y = vy[1] + m_v * (i_x - vx[1])
    return [round(i_x,const.n_round), round(i_y,const.n_round)]

#

def getIntersectionX(XX, zz):
    Xx = XX[0]
    #Xx = Xx[0]
    zx = zz[0]
    zy = zz[1]
    #zx1 = zx[0]
    #zx2 = zx[1]
    #zy1 = zy[0]
    #zy2 = zy[1]
    #m_z = ((zy[1] - zy[0]) / (zx[1] - zx[0]))
    i_x = Xx[0]
    i_y = zy[1] + ((zy[1] - zy[0]) / (zx[1] - zx[0])) * (i_x - zx[1])
    return [round(i_x,const.n_round), round(i_y,const.n_round)]

#
"""
def getAllIntersections(lines_):
    xlines = lines_[0]
    ulines = lines_[1]
    vlines = lines_[2]
    i_UV = getIntersectionUV(ulines, vlines)
    i_XU = getIntersectionX(xlines, ulines)
    i_XV = getIntersectionX(xlines, vlines)
    return [i_UV, i_XU, i_XV]
"""
############################################################## find coordinates
"""
def getMidpoint(lines_):
    allIntersections = getAllIntersections(lines_)
    i_UV = allIntersections[0]
    i_XU = allIntersections[1]
    i_XV = allIntersections[2]
    mp_x = (i_UV[0] + i_XU[0] + i_XV[0])/3
    mp_y = (i_UV[1] + i_XU[1] + i_XV[1])/3
    return [mp_x, mp_y]
"""
#

def getMidpoint(lines_):
    i_UV = getIntersectionUV(lines_[1], lines_[2])
    i_XU = getIntersectionX(lines_[0], lines_[1])
    i_XV = getIntersectionX(lines_[0], lines_[2])
    mp_x = (i_UV[0] + i_XU[0] + i_XV[0])/3
    mp_y = (i_UV[1] + i_XU[1] + i_XV[1])/3
    return [mp_x, mp_y]

#

def getMPorIntersection(wires):
    typematch = [1,1,1]
    if wires[0][0][0] == -999:
        typematch[2] = 0
    if wires[1][0][0] == -999:
        typematch[1] = 0
    if wires[2][0][0] == -999:
        typematch[0] = 0
    
    if sum(typematch) == 3:
        mp = getMidpoint(wires)
    elif sum(typematch) == 2:
        if typematch[2] == 0:
            mp = getIntersectionUV(wires[1], wires[2])
        elif typematch[1] == 0:
            mp = getIntersectionX(wires[0], wires[2])
        elif typematch[0] == 0:
            mp = getIntersectionX(wires[0], wires[1])
    else:
        mp = -666
        print("\n\nSOMETHING IS WRONG WITH MATCHING TRACKS???\n")
    return (mp, typematch)

#

def getBigTriangle(lines_):
    xlines = lines_[0]
    ulines = lines_[1]
    vlines = lines_[2]
    i_UV = getIntersectionUV(ulines, vlines)
    i1 = [i_UV[0], i_UV[1] + const.rr]
    i2 = [i_UV[0], i_UV[1] - const.rr]
    i_XU = getIntersectionX(xlines, ulines)
    i3 = [i_XU[0] + const.dd, i_XU[1] + const.rr/2]
    i4 = [i_XU[0] - const.dd, i_XU[1] - const.rr/2]
    i_XV = getIntersectionX(xlines, vlines)
    i5 = [i_XV[0] + const.dd, i_XV[1] - const.rr/2]
    i6 = [i_XV[0] - const.dd, i_XV[1] + const.rr/2]
    return [i1, i2, i3, i4, i5, i6]

#

def getHexagon(mp):
    x1 = mp[0] - const.dd
    x2 = mp[0]
    x3 = mp[0] + const.dd
    y1 = mp[1] + const.rr
    y2 = mp[1] + const.rr/2
    y3 = mp[1] - const.rr/2
    y4 = mp[1] - const.rr
    z1 = [x1, y3]
    z2 = [x1, y2]
    z3 = [x2, y1]
    z4 = [x3, y2]
    z5 = [x3, y3]
    z6 = [x2, y4]
    return [z1, z2, z3, z4, z5, z6]

#

def getWireCoords(wireplane, w_):
    # wire plane number = wireid + (338 - 332)/2 + 1
    # wireid 0 = wire plane number 4
    # wireid 331 = wire plane number 335
    if w_[wireplane] > 331 or (w_[wireplane] < 0 and w_[wireplane] != -999): print("uh oh: wire number confusion???")
    if w_[wireplane] == -999:
        x1 = -999
        x2 = -999
        y1 = -999
        y2 = -999
    else:
        wirenum = w_[wireplane] + 4
        wirenum = wirenum + 0.5
        #x plane
        if wireplane == 0:
            x = 2*const.dd * (wirenum - const.radius)
            x1 = x
            x2 = x
            y = const.radius
            y1 = - 1 * y
            y2 = y
        #u or v plane
        else:
            x = const.radius * 2*const.dd
            x1 = x
            x2 = -1 * x
            y1 = wirenum - 255
            y2 = wirenum - 85
            #v is sign difference from u
            if wireplane == 2:
                y1 = -1 * y1
                y2 = -1 * y2
    return [[x1, x2], [y1, y2]]

#

def isPointInPoly(nvert, vertx, verty, testx, testy):
    j = nvert-1
    c = False
    for i in range(nvert):
        if  ((verty[i]>testy) != (verty[j]>testy)) and (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) :
            c = ~c;
        j = i
    return c

################################################################### find values

def getTrackAngle(gmp_e,gmp_l):
    x_len = (gmp_l[0]-gmp_e[0])
    y_len = (gmp_l[1]-gmp_e[1])
    if x_len != 0:
        ang = math.atan(y_len/x_len)
    elif y_len > 0:
        ang = math.pi/2
    elif y_len < 0:
        ang = - 1 * math.pi/2    
    else:
        ang = -666
    return ang

#
"""
def PolygonArea(corners):
    #source: https://plot.ly/python/polygon-area/
    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area
"""
def PolygonArea(corners):
    ii = np.array(corners)
    jj = corners[1:len(corners)]
    jj.append(corners[0])
    jj = np.array(jj)
    corners = ii[:,0] * jj[:,1] - jj[:,0] * ii[:,1]
    area = abs(sum(corners))/2.0
    return round(area, const.n_round)

#

def getTrackLength(ints_e, ints_l):
    #thanks Eric Marzec
    coords_e = np.array(ints_e)
    coords_l = np.array(ints_l)
    x_e = coords_e[:,0]
    y_e = coords_e[:,1]
    x_l = coords_l[:,0]
    y_l = coords_l[:,1]
    x_matrix = x_e - x_l[:,np.newaxis]
    y_matrix = y_e - y_l[:, np.newaxis]
    dsqr = x_matrix**2 + y_matrix**2
    max_val = np.max(dsqr)
    return round(np.sqrt(max_val), const.n_round)

######################################################################### other

def combineAllPoints(ints_):
    # there is probably a MUCH better way to do this...
    ints_all = []
    for aa in range(len(ints_)):
        ints_all.append((ints_[aa][0], ints_[aa][1]))
    return ints_all

#

def convex_hull(points):
    ###########################################################################
    # Computes the convex hull of a set of 2D points.
    #
    # Input: an iterable sequence of (x, y) pairs representing the points.
    # Output: a list of vertices of the convex hull in counter-clockwise order,
    #   starting from the vertex with the lexicographically smallest coordinates.
    # Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    #
    # source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python
    ###########################################################################

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    return lower[:-1] + upper[:-1]

#

def enactLimitOnAllTracksInfo(alltracks, bal):
    goodtracks = allTrackInfo(alltracks.ts, alltracks.pf, alltracks.inf, bal, alltracks.ne)
    index = list(range(0,len(alltracks.ia_e)))
    ia_e = np.array(alltracks.ia_e)
    where = list(np.argwhere(ia_e > bal)[:,0])
    index = [aaa for aaa in index if aaa not in where]
    ia_l = np.array(alltracks.ia_l)    
    where = list(np.argwhere(ia_l > bal)[:,0])
    index = [aaa for aaa in index if aaa not in where]
    index = np.array(index)
    
    n_len = len(list(index))
    
    goodtracks.mt_e = list(np.array(alltracks.mt_e)[index])
    goodtracks.mt_l = list(np.array(alltracks.mt_l)[index])
    goodtracks.ia_e = list(np.array(alltracks.ia_e)[index])
    goodtracks.ia_l = list(np.array(alltracks.ia_l)[index])
    goodtracks.pa = list(np.array(alltracks.pa)[index])
    goodtracks.tl = list(np.array(alltracks.tl)[index])
    goodtracks.ang = list(np.array(alltracks.ang)[index])
    
    se = np.array(alltracks.se)[index]
    goodtracks.se = [[se[aaa][bbb] for bbb in range(len(se[aaa]))] for aaa in range(n_len)]
    w_e = np.array(alltracks.w_e)[index]
    goodtracks.w_e = [[w_e[aaa][bbb] for bbb in range(len(w_e[aaa]))] for aaa in range(n_len)]
    w_l = np.array(alltracks.w_l)[index]
    goodtracks.w_l = [[w_l[aaa][bbb] for bbb in range(len(w_l[aaa]))] for aaa in range(n_len)]
    mp_e = np.array(alltracks.mp_e)[index]
    goodtracks.mp_e = [[mp_e[aaa][bbb] for bbb in range(len(mp_e[aaa]))] for aaa in range(n_len)]
    mp_l = np.array(alltracks.mp_l)[index]
    goodtracks.mp_l = [[mp_l[aaa][bbb] for bbb in range(len(mp_l[aaa]))] for aaa in range(n_len)]
    ic_e = np.array(alltracks.ic_e)[index]
    
    goodtracks.ic_e = [[[ic_e[aaa][bbb][ccc] for ccc in range(len(ic_e[aaa][bbb]))] for bbb in range(len(ic_e[aaa]))] for aaa in range(n_len)]    
    ic_l = np.array(alltracks.ic_l)[index]
    goodtracks.ic_l = [[[ic_l[aaa][bbb][ccc] for ccc in range(len(ic_l[aaa][bbb]))] for bbb in range(len(ic_l[aaa]))] for aaa in range(n_len)]

    return goodtracks

###############################################################################
