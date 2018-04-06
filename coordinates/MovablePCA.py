#!/bin/python
import sys, os                            # environmental variables
sys.path.append('/home/jordancaraballo/Wolffia/')  # add wolffia repository for DCDReader utility
from lib.io.dcd import DCDReader          # include wolffia DCDReader library
import heapq                              # import for median
import math                               # import for distance calculation
import numpy as np                        # include numpy
import csv

#-------------------------------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------------------------------
mix      = "Tetradecane005"                           # name of the mixture - e.g Pentane09
workPath = "/media/jordancaraballo/b3cc799d-d1d9-43a9-a87d-321420b4a13b/home/Tetradecane/" + mix + "/"    # path to mixture - parental directory
psfFile  = workPath + mix + ".psf"                # name of the psf file
pdbFile  = workPath + mix + "_lastFrame.pdb"      # header of the pdb files
dcdFile  = workPath + mix + "_20ns.dcd" # dcd file

# Atom types to detect
upper_types = ["N0", "N1"] # atoms to get coordinates that we want
base_type   = 'O6'         # atom to calculate the height of pABA
mid_type    = 'C8'         # atom that connects pABA and Cellulose, pentane is C8, tetradecane is C15

# Split lines from coordinates, atom, and trajectories files
psf_list = [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]       # list with psf file lines
pdb_list = [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:]       # list with pdb file lines
dcd      = DCDReader(dcdFile) # list with dcd files

# Get coordinates from a single frame
def getCoordinatesSingleFrame(pdb_list, psf_list, singleDCD, molLength):
    tempCoors   = list() # list to store temporary coordinates for each pABA
    coordinates = list() # list to store final and clean coordinates
    for frame in singleDCD: # iterate over every frame
        for atom in range(len(pdb_list)): # iterate over every atom in the frame
            try:
                if psf_list[atom][3] == 'PAB': # if atom is PABA go for it
                    tempCoors.append(float(frame[0][int(pdb_list[atom][1])-1])) # append x coordinates
                    tempCoors.append(float(frame[1][int(pdb_list[atom][1])-1])) # append y coordinates
                    tempCoors.append(float(frame[2][int(pdb_list[atom][1])-1])) # append z coordinates
                    if len(tempCoors) == molLength: # if tempCoors has all the atoms from the molecule
                        coordinates.append(tempCoors) # send atoms to the main list
                        tempCoors = []                # empty temp list
            except IndexError: # some lines are trash lines from the pdb file
                pass
    return coordinates # multidimensional array

def getCoordinatesMultipleFrames(pdb_list, psf_list, multDCD, molLength, ranges):
    tempCoors   = list() # list to store temporary coordinates for each pABA
    coordinates = list() # list to store final and clean coordinates
    counter     = 0      # counter to increment frame periodically
    for frame in multDCD: # iterate over every frame
        #print counter, " % ", ranges, "  = ", counter % ranges
        if counter % ranges == 0:
            for atom in range(len(pdb_list)): # iterate over every atom in the frame
                try:
                    if psf_list[atom][3] == 'PAB': # if atom is PABA go for it
                        tempCoors.append(float(frame[0][int(pdb_list[atom][1])-1])) # append x coordinates
                        tempCoors.append(float(frame[1][int(pdb_list[atom][1])-1])) # append y coordinates
                        tempCoors.append(float(frame[2][int(pdb_list[atom][1])-1])) # append z coordinates
                        if len(tempCoors) == molLength: # if tempCoors has all the atoms from the molecule
                            coordinates.append(tempCoors) # send atoms to the main list
                            tempCoors = []                # empty temp list
                except IndexError: # some lines are trash lines from the pdb file
                    pass
        else:
            pass
        counter = counter + 1
    return coordinates # multidimensional array

def saveCoordinatesFile(coorsList, FileName):
    with open(FileName, 'wb') as myfile:
        wr = csv.writer(myfile)
        for i in coorsList:
            wr.writerow(i)

#cor = getCoordinatesSingleFrame(pdb_list, psf_list, dcd, 156)
coor = getCoordinatesMultipleFrames(pdb_list, psf_list, dcd, 156, 50)
saveCoordinatesFile(coor, "Tetradecane005MultDCDR50.txt")
print len(coor)
