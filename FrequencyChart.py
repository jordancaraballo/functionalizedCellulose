#!/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------------------------------
# Script developed to create a bubble chart from the frequencies and norms of files that hold
# coordinates. To run this script you will need files populated with numbers separated by
# end line separators. This code is written for a specific application, but with minor changes it
# can fit many needs.
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
#-------------------------------------------------------------------------------------------------------

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Global variables for the graph
name    = "TetradecaneHeights/Tetradecane"
files   = ["005", "01", "02","03","04","05","06","07","08","09","100"]
yValues = [5.00, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.00]
colors  = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','cyan', 'fuchsia', "orange"]

# Set array with data
# Filter and calculate frequencies
data    = []
minList = []
maxList = []

for i in range(len(files)):
    f = open(name + files[i] + ".txt") # open file with coordinates
    tmpArr = np.array([float(x) for x in f.read().split("\n")[:-1]]) # split and clean the data from the file
    print tmpArr
    minList.append(min(tmpArr))
    maxList.append(max(tmpArr))
    data.append(tmpArr)
dataParsed = np.asarray(data)


# delta values to produce the frequency
dMin = min(minList)
dMax = max(maxList)

# Define figures for the plot
fig = plt.figure()
fig.suptitle('Distribution of PABA Heights - Minus Base',fontsize=18)

ax  = fig.add_subplot(111)
fig.subplots_adjust(top=0.90)

ax.set_title('14-atoms spacer arm')

ax.set_xlabel('PABA Height (A)')
ax.set_ylabel('PABA Occupancy')

# Filter and calculate frequencies
for i in range(len(files)):

    # Calculate frequencies
    h = np.histogram(data[i],20,(dMin,dMax),density=True)

    n = len(h[1][:-1])
    x = h[1][:-1]
    y = np.ones(n)*yValues[i]
    sizes=(h[0]*5000)

    ax.scatter(x,y,s=sizes,c=colors[i])

plt.show()
