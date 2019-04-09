#!/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------------------------------
# Script developed to create a PCA chart from the frequencies and norms of files that hold
# coordinates. The coordinates are taken from frames at the desired range. Coordinates can be
# calculated using the MovablePCA.py script.
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
#-------------------------------------------------------------------------------------------------------
# Some background information
# Each Tetradecane molecule has 52 atoms
# Tetradecane005 has 156 dimensions in total in the x axis (3 molecules x 52 = 156)
# Tetradecane005 every 100 frames shape is (600, 156)
# Tetradecane01 every 100 frames shape is (1400, 156)
# Tetradecane02 every 100 frames shape is (2800, 156)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import datasets
from sklearn.decomposition import PCA
import numpy as np
import math

# Import tetradecane data
data005 = np.loadtxt('../data/TetradecanePCA/Tetradecane005MultDCDR100.txt',delimiter=',')
data01  = np.loadtxt('../data/TetradecanePCA/Tetradecane01MultDCDR100.txt',delimiter=',')
data02  = np.loadtxt('../data/TetradecanePCA/Tetradecane02MultDCDR100.txt',delimiter=',')
allData = np.concatenate((data005, data01, data02), axis=0)

print (data005.shape)
print (data01.shape)
print (data02.shape)
#print allData
#print " Shape " , allData.shape, allData.shape[0]
for i in range(allData.shape[0]):
    base = allData[i][:3].copy()
    for j in range(0, allData.shape[1], 3):
        allData[i,j] = allData[i,j] - base[0]
        allData[i,j+1] = allData[i,j+1] - base[1]
        allData[i,j+2] = allData[i,j+2] - base[2]
    #print allData[i,:10]

# To getter a better understanding of interaction of the dimensions
# plot the first three PCA dimensions
fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=110)
pca = PCA(n_components=3)
X_reduced =pca .fit_transform(allData)
#print "PCA=",pca.singular_values_

print (X_reduced.shape)
#print X_reduced[3:10, 0]
simMols = X_reduced[:600]

from operator import itemgetter

segregatedMols  = [simMols[0::3], simMols[1::3], simMols[2::3]]
distances = list()
rango = 100
for i in segregatedMols:
    tempDist = list()
    for index in range(len(i)-1):
        dist = math.sqrt((i[index][0] - i[index + 1][0])**2 + (i[index][1] - i[index + 1][1])**2 + (i[index][2] - i[index + 1][2])**2)
        tempDist.append([(index+1) * rango, dist])
    distances.append(sorted(tempDist, key=lambda x: x[1], reverse=True))

print (distances[0][:4])
print (distances[1][:4])
print (distances[2][:4])

# First Simulation
ax.scatter3D(simMols[0::3, 0], simMols[0::3, 1], simMols[0::3, 2], c='b')
ax.scatter3D(simMols[1::3, 0], simMols[1::3, 1], simMols[1::3, 2], c='g')
ax.scatter3D(simMols[2::3, 0], simMols[2::3, 1], simMols[2::3, 2], c='r')


ax.set_title("First three PCA directions")
ax.set_xlabel("1st eigenvector")
ax.w_xaxis.set_ticklabels([])
ax.set_ylabel("2nd eigenvector")
ax.w_yaxis.set_ticklabels([])
ax.set_zlabel("3rd eigenvector")
ax.w_zaxis.set_ticklabels([])

plt.show()

"""
Main code execution
"""
#if __name__ == "__main__":
