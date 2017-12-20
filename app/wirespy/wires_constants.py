# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:21:12 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: constants & geometry & some visual options for wires.py
"""

#import math

###############################################################################
# constants & geometry
###############################################################################
maxwires = 338 #this works for even numbers of wires only -- this describes the wire geometry

radius = maxwires/2 + 1

dd = 0.4330127019 #math.cos(math.pi/6)/2
rr = 0.5 #math.sin(math.pi/6)

xaxis = [[ -1 * radius * 2*dd, radius * 2*dd], [0, 0]]
uaxis = [[-1 * radius * dd, radius * dd], [-1 * radius * (rr + 1)/2, radius * (rr + 1)/2]]
vaxis = [[-1 * radius * dd, radius * dd], [radius * (rr + 1)/2, -1 * radius * (rr + 1)/2]]
#axis = [xaxis, uaxis, vaxis]

n_round = 10

###############################################################################
# visuals
###############################################################################
plotaxis = 0    #shows axis on hexagon

hotshading = 0.05   #alpha of shading
ms = 2
lw = 0.5

startcolor = (0, 255, 0)
startcolor = '#%02x%02x%02x' % startcolor
endcolor = (255, 0, 255)
endcolor = '#%02x%02x%02x' % endcolor
midcolor = (0, 128, 255)
midcolor = '#%02x%02x%02x' % midcolor        
