#!/bin/python
import sys, os                            # environmental variables
sys.path.append('/home/jordancaraballo/Wolffia/')  # add wolffia repository for DCDReader utility
from lib.io.dcd import DCDReader          # include wolffia DCDReader library
import heapq                              # import for median
import math                               # import for distance calculation
import numpy as np                        # include numpy
#-------------------------------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------------------------------
mix      = "Tetradecane100"                           # name of the mixture - e.g Pentane09
workPath = "/media/jordancaraballo/Datos_Frances/jordancaraballo/Investigacion/Tetradecane/" + mix + "/"    # path to mixture - parental directory
psfFile  = workPath + mix + ".psf"                # name of the psf file
pdbFile  = workPath + mix + "_lastFrame.pdb"      # header of the pdb files
dcdFile  = workPath + mix + "_LastFrame.dcd" # dcd file

# Atom types to detect
upper_types = ["N0", "N1"] # atoms to get coordinates that we want
base_type   = 'O6'         # atom to calculate the height of pABA
mid_type    = 'C8'         # atom that connects pABA and Cellulose, pentane is C8, tetradecane is C15

# Split lines from coordinates, atom, and trajectories files
psf_list = [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]       # list with psf file lines
pdb_list = [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:]       # list with pdb file lines
dcd      = DCDReader(dcdFile) # list with dcd files

z_coor = list()
for frame in dcd:
    for atom in range(len(pdb_list)): # for range in the amount of atoms
        try:
            # if atom is PABA go for it
            if psf_list[atom][3] == 'PAB':
                #print int(pdb_list[atom][1]) - 1
                z_coor.append(   float   (   frame[0][int(pdb_list[atom][1])-1]   )   )
                z_coor.append(   float   (   frame[1][int(pdb_list[atom][1])-1]   )   )
                z_coor.append(   float   (   frame[2][int(pdb_list[atom][1])-1]   )   )
        except IndexError:
            pass

#print z_coor
new = []
for i in range(0, len(z_coor), 156):
    new.append(z_coor[i : i+156])

for i in new:
    print i
"""
from sklearn.decomposition import PCA as sklearnPCA
sklearn_pca = sklearnPCA(n_components=2)
Y_sklearn = sklearn_pca.fit_transform(X_std)



traces = []

for name in ('Iris-setosa', 'Iris-versicolor', 'Iris-virginica'):

    trace = Scatter(
        x=Y_sklearn[y==name,0],
        y=Y_sklearn[y==name,1],
        mode='markers',
        name=name,
        marker=Marker(
            size=12,
            line=Line(
                color='rgba(217, 217, 217, 0.14)',
                width=0.5),
            opacity=0.8))
    traces.append(trace)


data = Data(traces)
layout = Layout(xaxis=XAxis(title='PC1', showline=False),
                yaxis=YAxis(title='PC2', showline=False))
fig = Figure(data=data, layout=layout)
py.iplot(fig)
"""
