# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:21:12 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: constants & geometry & some visual options for wires.py
"""

import math

###############################################################################
# constants & geometry
###############################################################################
maxwires = 338 #this works for even numbers of wires only -- this describes the wire geometry

radius = maxwires/2 + 1

xaxis = [[ -1 * radius * math.cos(math.pi/6), radius * math.cos(math.pi/6)], [0, 0]]
uaxis = [[-1 * radius * math.cos(math.pi/6)/2, radius * math.cos(math.pi/6)/2], [-1 * radius * (math.sin(math.pi/6) + 1)/2, radius * (math.sin(math.pi/6) + 1)/2]]
vaxis = [[-1 * radius * math.cos(math.pi/6)/2, radius * math.cos(math.pi/6)/2], [radius * (math.sin(math.pi/6) + 1)/2, -1 * radius * (math.sin(math.pi/6) + 1)/2]]
#axis = [xaxis, uaxis, vaxis]


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
