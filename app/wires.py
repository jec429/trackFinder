# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:21:32 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description:
copied & working off of wire.py created by Jorge Chaves
newest change: slimmed down unused things + prepping for fixed coordinates
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
import os,sys

import ROOT

###############################################################################
# features to toggle
###############################################################################
plotaxis = 0    #shows axis on hexagon
plotEach = 0

mapintersection = 0 #0=no shading, 1=shade area between 3 wires
hotshading = 0.05   #alpha of shading
badarealimit = 5000

###############################################################################
# constants & geometry
###############################################################################
maxwires = 338 #this works for even numbers of wires only -- this describes the wire geometry

radius = maxwires/2 + 1
c_pi6 = math.cos(math.pi/6)
s_pi6 = math.sin(math.pi/6)

xaxis = [[ -1 * radius * c_pi6, radius * c_pi6], [0, 0]]
uaxis = [[-1 * radius * c_pi6/2, radius * c_pi6/2], [-1 * radius * (s_pi6 + 1)/2, radius * (s_pi6 + 1)/2]]
vaxis = [[-1 * radius * c_pi6/2, radius * c_pi6/2], [radius * (s_pi6 + 1)/2, -1 * radius * (s_pi6 + 1)/2]]
axis = [xaxis, uaxis, vaxis]

startcolor = (0, 255, 0)
startcolor = '#%02x%02x%02x' % startcolor
endcolor = (255, 0, 255)
endcolor = '#%02x%02x%02x' % endcolor
midcolor = (0, 128, 255)
midcolor = '#%02x%02x%02x' % midcolor        



###############################################################################
# functions
###############################################################################

############################################################ plotting funtcions

def plotAxis():
    color = (200, 200, 200)
    color = '#%02x%02x%02x' % color
    plt.plot(xaxis[0], xaxis[1], color, ls='-')
    plt.plot(uaxis[0], uaxis[1], color, ls='-')
    plt.plot(vaxis[0], vaxis[1], color, ls='-')

#

def plotMPs(mp_e, mp_l):
    plt.plot(mp_e[0], mp_e[1], startcolor, marker='.')
    plt.plot(mp_l[0], mp_l[1], endcolor, marker='.')

#

def plotAllRegions(allmps, allints_all, prefix):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, aspect='equal')
    ax1.add_patch(
        patches.RegularPolygon(
            (0, 0),      # (x,y)
            6,           # number of vertices
            radius,      # radius
            fill = False
        )
    )
    ax1.set_xlim(-200,200)
    ax1.set_ylim(-200,200)
    if plotaxis == 1: plotAxis()
    allmps_e = allmps[0]
    allmps_l = allmps[1]
    for aa in range(len(allmps_e)):
        mp_e = allmps_e[aa]
        mp_l = allmps_l[aa]
        ints_all = allints_all[aa]
        
        if mapintersection == 1:
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_all, color=midcolor, alpha=hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
        
        plt.plot([mp_e[0], mp_l[0]], [mp_e[1], mp_l[1]], midcolor, ls='-')
        plotMPs(mp_e, mp_l)
    plt.savefig('plots/'+sys.argv[1].replace('.root','')+'/'+prefix+'midpoints.pdf', bbox_inches='tight')
    plt.savefig('plots/'+sys.argv[1].replace('.root','')+'/'+prefix+'midpoints.png', bbox_inches='tight')

################################################################### find values

def getWireCoords(wireplane, w_):
    # wire plane number = wireid + (338 - 332)/2 + 1
    # wireid 0 = wire plane number 4
    # wireid 331 = wire plane number 335
    if w_[wireplane] > 331: print("uh oh: wire number confusion???")
    wirenum = w_[wireplane] + 4
    wirenum = wirenum + 0.5
    #x plane
    if wireplane == 0:
        x = c_pi6 * (wirenum - radius)
        x1 = x
        x2 = x
        y = radius
        y1 = - 1 * y
        y2 = y
    #u or v plane
    else:
        x = radius * c_pi6
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

def getMidpoint(lines_):
    allIntersections = getAllIntersections(lines_)
    i_UV = allIntersections[0]
    i_XU = allIntersections[1]
    i_XV = allIntersections[2]
    mp_x = (i_UV[0] + i_XU[0] + i_XV[0])/3
    mp_y = (i_UV[1] + i_XU[1] + i_XV[1])/3
    r_error = math.sqrt((i_UV[0] - mp_x)**2 + (i_UV[1] - mp_y)**2)
    return [mp_x, mp_y, r_error]

############################################################ find intersections

def getAxisIntersection(wireplane, xx, yy):
    x1 = xx[0]
    x2 = xx[1]
    y1 = yy[0]
    y2 = yy[1]
    ax = axis[wireplane]
    if wireplane == 0:
        i_x = x1
        i_y = 0
    else:
        slope = (y2 - y1) / (x2 - x1)
        i_x = (y2 - slope * x2) / (ax[1][0]/ax[0][0] - slope)
        i_y = ax[1][0]/ax[0][0] * i_x
    dist = math.sqrt(i_x * i_x + i_y * i_y)
    return [i_x, i_y, dist]

#

def getIntersectionUV(uu, vv):
    ux = uu[0]
    uy = uu[1]
    vx = vv[0]
    vy = vv[1]
    ux1 = ux[0]
    ux2 = ux[1]
    uy1 = uy[0]
    uy2 = uy[1]
    m_u = (uy2 - uy1) / (ux2 - ux1)
    vx1 = vx[0]
    vx2 = vx[1]
    vy1 = vy[0]
    vy2 = vy[1]
    m_v = (vy2 - vy1) / (vx2 - vx1)
    i_x = (uy2 - m_u * ux2 - vy2 + m_v * vx2) / (m_v - m_u)
    i_y = vy2 + m_v * (i_x - vx2)
    return [i_x, i_y]

#

def getIntersectionX(XX, zz):
    Xx = XX[0]
    Xx = Xx[0]
    zx = zz[0]
    zy = zz[1]
    zx1 = zx[0]
    zx2 = zx[1]
    zy1 = zy[0]
    zy2 = zy[1]
    m_z = (zy2 - zy1) / (zx2 - zx1)
    i_x = Xx
    i_y = zy2 + m_z * (i_x - zx2)
    return [i_x, i_y]

#

def getAllIntersections(lines_):
    xlines = lines_[0]
    ulines = lines_[1]
    vlines = lines_[2]
    i_UV = getIntersectionUV(ulines, vlines)
    i_XU = getIntersectionX(xlines, ulines)
    i_XV = getIntersectionX(xlines, vlines)
    return [i_UV, i_XU, i_XV]

######################################################## unsorted new functions

def combineAllPoints(ints_e, ints_l):
    ints_all = []
    for aa in range(3):
        p1 = ints_e[aa]
        p1_x = p1[0]
        p1_y = p1[1]
        pp1 = (p1_x, p1_y)
        ints_all.append(pp1)
        p2 = ints_l[aa]
        p2_x = p2[0]
        p2_y = p2[1]
        pp2 = (p2_x, p2_y)
        ints_all.append(pp2)
    return ints_all

#

def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    
    source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#Python
    """

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

def PolygonArea(corners):
    """source: https://plot.ly/python/polygon-area/
    """
    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

#



###############################################################################
# main function
###############################################################################

def wiresMain(ew,lw,subevent,prefix):
    allmps_e = []
    allmps_l = []
    allints_all = []
    allareas = []
    badevents = []
    progress = 0

    h_angle = ROOT.TH1F("h_angle","",100,-3.15,3.15)
    
    for w_e,w_l,se in zip(ew,lw,subevent):
        lines_e = []
        lines_l = []
        
        if plotEach == 1:
            #plot polygon
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(
                patches.RegularPolygon(
                    (0, 0),      # (x,y)
                    6,           # number of vertices
                    radius,      # radius
                    fill = False
                )
            )
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
            if plotaxis == 1: plotAxis()
        
        #plot 3 wire planes
        for aa in range(3):
            coords = getWireCoords(aa, w_e)
            xx = coords[0]
            yy = coords[1]
            if plotEach == 1:
                plt.plot(xx, yy, startcolor, ls='-')
                        
            lines_e.append([xx,yy])
            
            coords = getWireCoords(aa, w_l)
            xx = coords[0]
            yy = coords[1]
            if plotEach == 1:
                plt.plot(xx, yy, endcolor, ls='-')

            lines_l.append([xx,yy])
        
        mp_e = getMidpoint(lines_e)
        mp_l = getMidpoint(lines_l)
             
        ints_e = getAllIntersections(lines_e)
        ints_l = getAllIntersections(lines_l)

        ints_all = combineAllPoints(ints_e, ints_l)
        ints_all = convex_hull(ints_all)
        ints_all = np.array(ints_all)
        
        allareas.append(PolygonArea(ints_e))
        allareas.append(PolygonArea(ints_l))
        
        if PolygonArea(ints_e) > badarealimit or PolygonArea(ints_l) > badarealimit:
            badevents.append(se)
        else:
            allmps_e.append(mp_e)
            allmps_l.append(mp_l)
            allints_all.append(ints_all)

        if plotEach == 1:
            if mapintersection == 1:
                ax1 = fig1.add_subplot(111, aspect='equal')
                ax1.add_patch(patches.Polygon(ints_all, color=midcolor, alpha=hotshading, lw=0))
                ax1.set_xlim(-200,200)
                ax1.set_ylim(-200,200)
                
                plotMPs(mp_e, mp_l)
            
            plt.savefig('plots/events/'+sys.argv[1].replace('.root','')+'/'+prefix+str(se)+'.pdf', bbox_inches='tight')
            plt.savefig('plots/events/'+sys.argv[1].replace('.root','')+'/'+prefix+str(se)+'.png', bbox_inches='tight')
            if progress % 25 == 0: print(prefix + str(se) + ".pdf saved\t~ " + str(100*progress/len(subevent)) + "% done")
            plt.close()
        else:
            if progress % 25 == 0: print("\t~ " + str(100*progress/len(subevent)) + "% done")

        progress = progress + 1
        
    plotAllRegions([allmps_e, allmps_l], allints_all, prefix)
    for gmp_e,gmp_l in zip(allmps_e, allmps_l):
        x_len = (gmp_l[0]-gmp_e[0])
        y_len = (gmp_l[1]-gmp_e[1])
        h_angle.Fill(math.atan(y_len/x_len))
    
    c1 = ROOT.TCanvas("c1","",900,600)
    h_angle.Draw()
    c1.Print('plots/'+sys.argv[1].replace('.root','')+'/'+'h_angle_'+prefix[:-1]+'.pdf')
    c1.Print('plots/'+sys.argv[1].replace('.root','')+'/'+'h_angle_'+prefix[:-1]+'.png')
    
#    print("\n\n\n\n\nbad events = ")
#    print badevents
#    print("\n\n\n\n\n")
#    return allints_all
    return allareas


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


# Read from ROOT file
if len(sys.argv) < 2:
    raise ValueError('Missing input file.')
    
f = ROOT.TFile(sys.argv[1])
t = ROOT.TTree()
f.GetObject("tracks",t)
os.system('mkdir -p plots/'+sys.argv[1].replace('.root',''))

ev = 1
for e in t:
    for i in xrange(0,e.nmatches):
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


"""
print("\n\new_beam = ")
print ew_beam
print("\n\nlw_beam = ")
print lw_beam
print("\n\nsubevent_beam = ")
print subevent_beam
print("\n\new_cosmics = ")
print ew_cosmics
print("\n\nlw_cosmics = ")
print lw_cosmics
print("\n\nsubevent_cosmics = ")
print subevent_cosmics
print("\n\n\n")
"""

###############################################################################
# run
###############################################################################

ans1 = wiresMain(ew_beam,lw_beam,subevent_beam,'plot_beam_wires_')
ans2 = wiresMain(ew_cosmics,lw_cosmics,subevent_cosmics,'plot_cosmics_wires_')




