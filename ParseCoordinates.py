#!/bin/python
#-------------------------------------------------------------------------------------------------------
# Script developed to extract the coordinates from certain type of atoms from simulations. This script needs
# the .psf file and pdb files that will store the coordinates (refer to the from_dcd_to_pdb.sh script).
# This script can be used to get coordinates from specific atoms in a trendign way from where the user
# can calculate height, frequency, histograms, PCA's, and many other analysis. To run this script the user
# will need to have the path to the pdb files and a psf from the same simulation, together with the type of atoms.
# For now, you can run the script as python ./PCA_CelluloseFrames.py > nameFile.txt
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
#-------------------------------------------------------------------------------------------------------
import heapq
import numpy as np

# Global variables
workPath  = "None"
psfFile   = workPath + "name.psf"  # name of the psf file - you can also add the path like "path/name.psf"
headFile  = "Example"              # header of the pdb files - will work perfectly if you used from_dcd_to_pdb script
InitFrame = 9001                   # first frame or number of the file - it is later on appended to the name
LastFrame = 10000                  # last frame or number of the file - will determine the lenght of the loop

# Atom types to detect
ring_types  = ["C0", "C1", "C2", "C3", "C4"]
upper_types = ["N0", "N1"]

# Creates list of every pdb file that will be located in the same path
coor_files = list()
for i in range(InitFrame,LastFrame+1,1):
    coor_files.append(headFile + str(i) + ".pdb")

# Creates list of each line on the psf file - this includes each atom info
psf_list = [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]

# Multidemensional list with the coordinates of the files
coor_list = list()
for frame in coor_files:
    tempList = list()
    for line in open(workPath + frame,'r').readlines():
        tempList.append(line.strip().split())
    coor_list.append(tempList[1:])

#-------------------------------------------------------------------------------------
# Search for PABA two top atoms (two external Nitrogens) and get a mean of them.
# Then do a dendogram of each of them that would express their height.
# The coor to get is in the z axis.
#--------------------------------------------------------------------------------------
heightCoor_list = list()
z_coor          = list()

for frame in coor_list:
    z_coor_list = list()
    cont = 0
    for atom in range(len(frame)):
        try:
            # if this atom is part of the selected ones
            if (psf_list[atom][5])[:2] in upper_types:
                z_coor.append(float(frame[atom][8])) # append to list
                cont += 1 # add to the counter
            # if they passed the amount of Nitrogen atoms then save it or print it
            if cont == 3:
                print np.array(heapq.nlargest(2,z_coor)).mean() # print mean to store it in a file
                z_coor = [] # empty list
                cont = 0    # renew counter

        except IndexError:
            pass
