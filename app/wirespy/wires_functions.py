# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:12:36 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: functions that wires.py uses to plot
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import math

import wires_constants as const
import wires_toggle as tog



############################################################ plotting functions

def plotAxis():
    color = (200, 200, 200)
    color = '#%02x%02x%02x' % color
    plt.plot(const.xaxis[0], const.xaxis[1], color, ls='-')
    plt.plot(const.uaxis[0], const.uaxis[1], color, ls='-')
    plt.plot(const.vaxis[0], const.vaxis[1], color, ls='-')

#

def plotMPs(mp_e, mp_l):
    plt.plot(mp_e[0], mp_e[1], const.startcolor, marker='.')
    plt.plot(mp_l[0], mp_l[1], const.endcolor, marker='.')

#

def plotAllRegions(allmps, polyregion_all, prefix):
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
    for aa in range(len(allmps_e)):
        mp_e = allmps_e[aa]
        mp_l = allmps_l[aa]
        ints_all = polyregion_all[aa]
        
        if tog.plt_shading == 1:
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_all, color=const.midcolor, alpha=const.hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
        
        if tog.plt_lines == 1: plt.plot([mp_e[0], mp_l[0]], [mp_e[1], mp_l[1]], const.midcolor, ls='-')
        if tog.plt_dots == 1: plotMPs(mp_e, mp_l)
    if tog.locallaptop == 1:
        plt.savefig(prefix+'midpoints.pdf', bbox_inches='tight')
    else:
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
        x = math.cos(math.pi/6) * (wirenum - const.radius)
        x1 = x
        x2 = x
        y = const.radius
        y1 = - 1 * y
        y2 = y
    #u or v plane
    else:
        x = const.radius * math.cos(math.pi/6)
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
    
############################################################ find intersections

def getAxisIntersection(wireplane, xx, yy):
    x1 = xx[0]
    x2 = xx[1]
    y1 = yy[0]
    y2 = yy[1]
    ax = const.axis[wireplane]
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

######################################################################### other

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

###################################################################new/unsorted
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
        print("\n\nSOMETHING VERY WEIRD IS GOING ON WITH ANGLES!")
        print("x_len = " + str(x_len) + "\ty_len = "+ str(y_len))
        ang = -666
    return ang

#