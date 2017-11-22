# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:08:34 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description:
copied & working off of wire.py created by Jorge Chaves
newest change: fixed wire coords yet again
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

###############################################################################
# features to toggle
###############################################################################
plotaxis = 1    #shows axis on hexagon
printstuff = 0  #prints out some info I was curious about at times while writing

mapintersection = 1 #0=no shading, 1=shade area between 3 wires
hotshading = 0.1   #alpha of shading

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

def getColor(time, t_e):
    t_min = min(t_e)
    t_max = t_min + 375
    t = (time - t_min)/(t_max - t_min)
    n_colors = 4
    if t > 2/(n_colors - 1):
        c_t = int((n_colors - 1) * (t - 2/(n_colors - 1)) * 255)
        color = (c_t, 0, 255)
    elif t > 1/(n_colors - 1):
        c_t = int((n_colors - 1) * (t - 1/(n_colors - 1)) * 255)
        color = (0, 255 - c_t, 255)
    else:
        c_t = int((n_colors - 1) * t * 255)
        color = (0, 255, c_t)
    color = '#%02x%02x%02x' % color
    return color

############################################################ plotting funtcions
"""
#currently unused
def plotPolygon():
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
"""
#

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

def plotAllMPs(allmps, allints):
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
    allints_e = allints[0]
    allints_l = allints[1]
    for aa in range(len(subevent)):
        mp_e = allmps_e[aa]
        mp_l = allmps_l[aa]
        ints_e = allints_e[aa]
        ints_l = allints_l[aa]
        
        if mapintersection == 1:
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_e, color=startcolor, alpha=hotshading, lw=0))
            ax1.add_patch(patches.Polygon(ints_l, color=endcolor, alpha=hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
        
        plt.plot([mp_e[0], mp_l[0]], [mp_e[1], mp_l[1]], midcolor, ls='-')
        plotMPs(mp_e, mp_l)
    plt.savefig('midpoints.pdf', bbox_inches='tight')


#
"""
#unused
def plotIntersections(xx, uu, vv):
    i_ = getIntersectionUV(uu, vv)
    i_x = i_[0]
    i_y = i_[1]
    plt.plot(i_x, i_y, 'black', marker='.')
    
    i_ = getIntersectionX(xx, uu)
    i_x = i_[0]
    i_y = i_[1]
    plt.plot(i_x, i_y, 'black', marker='.')

    i_ = getIntersectionX(xx, vv)
    i_x = i_[0]
    i_y = i_[1]
    plt.plot(i_x, i_y, 'black', marker='.')
"""
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
        y1 = wirenum - 85
        y2 = wirenum - 255
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

################################################################### prints outs

def printOutInfo(wireplane, w_, xx, yy):
    wirenum = w_[wireplane]
    if wireplane == 0:
        plane = "X"
    elif wireplane == 1:
        plane = "U"
    elif wireplane == 2:
        plane = "V"
    x1 = xx[0]
    x2 = xx[1]
    y1 = yy[0]
    y2 = yy[1]
    intersect = getAxisIntersection(wireplane, xx, yy)
    i_x = intersect[0]
    i_y = intersect[1]
    dist = intersect[2]
#    plt.plot(i_x, i_y, 'r*')
    print(plane + "\t" + str(wirenum) + "\t" + str(round(x1,2)) + "\t" + str(round(x2,2)) + "\t" + str(round(y1,2)) + "\t" + str(round(y2,2)) + "\t" + str(round(i_x,2)) + "\t" + str(round(i_y,2)) + "\t" + str(round(dist,2)))

#

def doChecks(t_e, t_l, se):
    #how color code?
    t_min = min(t_e)
    t_max = max(t_l)
    if t_max - t_min > 375:
        print(str(t_max - t_min))

    if max(t_e) >= min(t_l):
        print("BAD TIME: max(t_e) >= min(t_l)\t" + str(max(t_e)) + " >= " + str(min(t_l)) + "\tse = " + str(se))
    if t_l[0] < t_e[0]:
        print("BAD TIME 0!\t" + str(t_l[0]) + " < " + str(t_e[0]) + "\tse = " + str(se))
    if t_l[1] < t_e[1]:
        print("BAD TIME 1!\t" + str(t_l[1]) + " < " + str(t_e[1]) + "\tse = " + str(se))
    if t_l[2] < t_e[2]:
        print("BAD TIME 2!\t" + str(t_l[2]) + " < " + str(t_e[2]) + "\tse = " + str(se))

#






###############################################################################
# main function
###############################################################################

def wiresMain(ew,et,lw,lt,subevent,prefix):
    allmps_e = []
    allmps_l = []
    allints_e = []
    allints_l = []
    for w_e,t_e,w_l,t_l,se in zip(ew,et,lw,lt,subevent):
        lines_e = []
        lines_l = []
        
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

        if printstuff == 1: print("PLANE\twirenum\tx1\tx2\ty1\ty2\ti_x\ti_y\tdist")
        
        #plot 3 wire planes
        for aa in range(3):
            color = getColor(t_e[aa], t_e)
            coords = getWireCoords(aa, w_e)
            xx = coords[0]
            yy = coords[1]
            plt.plot(xx, yy, color, ls='--')
                        
            lines_e.append([xx,yy])

            if printstuff == 1: printOutInfo(aa, w_e, xx, yy)
            
            color = getColor(t_l[aa], t_e)
            coords = getWireCoords(aa, w_l)
            xx = coords[0]
            yy = coords[1]
            plt.plot(xx, yy, color, ls='-')

            lines_l.append([xx,yy])

            if printstuff == 1: printOutInfo(aa, w_l, xx, yy)
                
        # checks
        doChecks(t_e, t_l, se)
        
        mp_e = getMidpoint(lines_e)
        ints_e = getAllIntersections(lines_e)
        allmps_e.append(mp_e)
        allints_e.append(ints_e)
        mp_l = getMidpoint(lines_l)
        ints_l = getAllIntersections(lines_l)
        allmps_l.append(mp_l)
        allints_l.append(ints_l)
        
        if mapintersection == 1:
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(patches.Polygon(ints_e, color=startcolor, alpha=hotshading, lw=0))
            ax1.add_patch(patches.Polygon(ints_l, color=endcolor, alpha=hotshading, lw=0))
            ax1.set_xlim(-200,200)
            ax1.set_ylim(-200,200)
            
            plotMPs(mp_e, mp_l)
        
        plt.savefig(prefix+str(se)+'.pdf', bbox_inches='tight')
        
    plotAllMPs([allmps_e, allmps_l], [allints_e, allints_l])
    return [allmps_e, allmps_l]



###############################################################################
# data
###############################################################################

# [x,u,v]
et = []
ew = []

lt = []
lw = []

subevent = []

# check
"""
et = [[0,150,300], [0,150,300]]
ew = [[165,165,165], [0,0,0]]

lt = [[75,225,375], [75,225,375]]
lw = [[166,166,166], [331,331,331]]

subevent = [0.0,0.1]

#

et.append([0,150,300])
ew.append([165,165,165])
lt.append([75,225,375])
lw.append([166,166,331])
subevent.append(0.2)

#

et.append([375,187,0])
ew.append([0,0,0])
lt.append([375,187,0])
lw.append([331,331,331])
subevent.append(0.3)

et.append([375,187,0])
ew.append([165,165,165])
lt.append([375,187,0])
lw.append([166,166,166])
subevent.append(0.4)
"""
# real data
"""
et.append([813.213,826.996,827.613])
ew.append([168,    156,    150])
lt.append([1128.39,1094.76,1143.94])
lw.append([57,     160,    198])
subevent.append(5)

#

et.append([-2218.03,-2258.38,-2219.13])
ew.append([198,     92,      203])
lt.append([-1897.63,-2095.31,-1911.08])
lw.append([316,     53,      169.5])
subevent.append(6.0)

et.append([604.196,588.31,600.826])
ew.append([78,     181,   241])
lt.append([617.196,610.869,648.221])
lw.append([49,     198,    248])
subevent.append(6.1)

et.append([4580.94,4588.52,4577.6])
ew.append([129,    103,    283])
lt.append([4733.44,4683.01,4728.77])
lw.append([11,     198,    259])
subevent.append(6.2)

#

et.append([-2586.3, -2615.83,-2590.63])
ew.append([273,     34,      209])
lt.append([-2360.57,-2452.99,-2364.91])
lw.append([294,     46,      125])
subevent.append(7.0)

et.append([-640.471,-648.045,-641.663])
ew.append([200,     1,       289])
lt.append([-400.792,-507.73, -497.268])
lw.append([213,     29,      254])
subevent.append(7.1)

et.append([262.468,255.302,262.224])
ew.append([225,    134,    134])
lt.append([308.328,300.672,301.102])
lw.append([267,    129,    102])
subevent.append(7.2)

#

et.append([-936.465,-932.741,-909.501])
ew.append([95,      198,     208])
lt.append([-675.623,-747.147,-686.214])
lw.append([100,     188,     203])
subevent.append(9)

#

et.append([109.59,195.499,107.096])
ew.append([318,   95,     67])
lt.append([224.005,217.991,219.176])
lw.append([303,    80,     102])
subevent.append(11)

#

et.append([-2621.97,-2639.12,-2619.02])
ew.append([315,     302,     283])
lt.append([-2476.02,-2467.49,-2443.65])
lw.append([315,     214,     293])
subevent.append(12.0)

et.append([5410.18,5409.2, 5415.95])
ew.append([239,    175,    15])
lt.append([5475.58,5520.72,5484.94])
lw.append([251,    187,    9])
subevent.append(12.1)

et.append([5440.38,5465.67,5462.97])
ew.append([145,    230,    216])
lt.append([5472.57,5473.4, 5485.62])
lw.append([133,    241,    212])
subevent.append(12.2)
"""
#
"""
et.append([])
ew.append([])
lt.append([])
lw.append([])
subevent.append()
"""


# Jorge's listed subevents (with made-up times)
ew = [[184,156,150], [78,181,241], [224,134,134], [318,95,67]]
et = [[0,0,0],       [0,0,0],      [0,0,0],       [0,0,0]]

lw = [[104,178,216], [49,198,255], [267,129,102], [303,80,102]]
lt = [[375,375,375], [375,375,375],[375,375,375], [375,375,375]]

subevent = [5,6,7,11]
"""
# my attempted matches of Jorge's listed subevents (slightly different?)
et.append([813.213,826.996,827.613])
ew.append([168,    156,    150])
lt.append([1128.39,1094.76,1143.94])
lw.append([57,     160,    198])
subevent.append(5)

et.append([604.196,588.31,600.826])
ew.append([78,     181,   241])
lt.append([617.196,610.869,648.221])
lw.append([49,     198,    248])
subevent.append(6.1)

et.append([262.468,255.302,262.224])
ew.append([225,    134,    134])
lt.append([308.328,300.672,301.102])
lw.append([267,    129,    102])
subevent.append(7.2)

et.append([109.59,195.499,107.096])
ew.append([318,   95,     67])
lt.append([224.005,217.991,219.176])
lw.append([303,    80,     102])
subevent.append(11)
"""

###############################################################################
# run
###############################################################################

ans = wiresMain(ew,et,lw,lt,subevent,'12096_')


