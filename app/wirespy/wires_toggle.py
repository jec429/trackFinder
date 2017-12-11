# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:55:45 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: features of wires.py to toggle
"""

locallaptop = 0 #1=shannon's working on local laptop, 0=not that

plotEach = 0    #1=plot each individual track, 0=only plot visual heat maps 
badarealimit = 2000 #38841 #limit of triangular area for track to be included

#plot style options
plt_dots = 1
plt_lines = 0
plt_shading = 1

#return options for wiresMain
returnoption = 2 #0=polyregion_all, 1=triareas_all, 2=badevents

showprogress = 0 #1=print out percent done with wiresMain, 0=don't

createlog = 1 #1=create a text file with info, 0=don't
