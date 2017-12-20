# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:55:45 2017

@author: Shannon Glavin - sglavin@sas.upenn.edu (Python 3.4)

Description: features of wires.py to toggle
"""

locallaptop = 0 #1=shannon's working on local laptop, 0=not that

plotEach = 0    #1=plot each individual track, 0=only plot visual heat maps 

#plot style options
plt_lines = 1
plt_shading = 0 #a lot slower than lines

plt_angles = 0
saveplotformat = 0 #0=png, 1=pdf, 2=both


# currently being phased out of code:
badarealimit = -1 #limit of triangular area for track to be included, -1=no limit
createlog = 1 #1=create a text file with info, 0=don't
plt_dots = 1
