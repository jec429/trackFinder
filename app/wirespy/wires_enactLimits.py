# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 13:02:34 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description:
"""

import os,sys

import wires_toggle as tog
import wires_functions as func

#from wires_functions import allTrackInfo



###############################################################################
# initial start/setup
###############################################################################

if tog.locallaptop == 0:
    inputfile = sys.argv[1]
    limit = eval(sys.argv[2])
else:
    inputfile = "mumbojumbo/test.root"
    limit = 0.434
inputfile = inputfile.split('/')[-1]
os.system('mkdir -p plots/'+inputfile.replace('.root','')+'/limit'+str(limit).replace('.','_')+'/')

###############################################################################

def enactLimit(inf, pf, bal):
    alltracks = func.readAllTrackInfo(inf, pf)
    goodtracks = func.enactLimitOnAllTracksInfo(alltracks, bal)
    
    func.writeAllTrackInfo_1(goodtracks)
    func.writeAllTrackInfo_3(goodtracks, len(alltracks.se))

    func.writeTrackInfo(goodtracks)
    if tog.plt_lines == 1 or tog.plt_shading == 1:
        func.plotAllRegions(goodtracks)
    
    return goodtracks

###############################################################################
# run
###############################################################################

tracks_beam = enactLimit(inputfile,'beam',limit)
tracks_cosmics = enactLimit(inputfile,'cosmics',limit)
