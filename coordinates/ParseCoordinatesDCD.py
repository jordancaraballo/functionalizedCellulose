#!/bin/python
"""
# Script developed to extract the coordinates from certain type of atoms from simulations. This script needs
# the .psf, .pdb and a .dcd file that will store the coordinates (refer to the from_dcd_to_pdb.sh script).
# This script can be used to get coordinates from specific atoms in a trending way from where the user
# can calculate height, frequency, histograms, PCA's, and many other analysis. To run this script the user
# will need to have the path to the pdb and psf files from the same simulation, together with the type of atoms.
# For now, you can run the script as python ./ParseCoordinates.py nameFile.dcd > nameFile.txt
# To create the pdb or dcd files you just need to run:
#
#   catdcd -o Tetradecane100_Last1000Frames.dcd -first 9001 -last 10000 Tetradecane100.dcd
#   catdcd -o Tetradecane08_lastFrame.pdb -otype pdb -s Tetradecane08.pdb -stype pdb -first 10000 -last 10000 Tetradecane08.dcd
#
# This program needs a library from the wolffia project that can be cloned from https://github.com/compMathUPRH/wolffia
#
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
"""
#-------------------------------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------------------------------
import sys, os                            # environmental variables
sys.path.append('/home/jordancaraballo/Wolffia/')  # add wolffia repository for DCDReader utility
from lib.io.dcd import DCDReader          # include wolffia DCDReader library
import heapq                              # import for median
import math                               # import for distance calculation
import numpy as np                        # include numpy
#-------------------------------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------------------------------
mix      = "Tetradecane02"                           # name of the mixture - e.g Pentane09
workPath = "/media/jordancaraballo/Datos_Frances/jordancaraballo/Investigacion/Tetradecane/" + mix + "/"    # path to mixture - parental directory
psfFile  = workPath + mix + ".psf"                # name of the psf file
pdbFile  = workPath + mix + "_lastFrame.pdb"      # header of the pdb files
dcdFile  = workPath + mix + "_Last1000Frames.dcd" # dcd file

# Atom types to detect
upper_types = ["N0", "N1"] # atoms to get coordinates that we want
base_type   = 'O6'         # atom to calculate the height of pABA
mid_type    = 'C8'         # atom that connects pABA and Cellulose, pentane is C8, tetradecane is C15

# Split lines from coordinates, atom, and trajectories files
psf_list = [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]       # list with psf file lines 
pdb_list = [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:]       # list with pdb file lines
dcd = [DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile)] # list with dcd files
#------------------------------------------------------------------------
"""
 ###################################################################################
 Function: getIndexUppAndMidTypes(frames, pdb_list, psf_list, upper_types, mid_type)
 ###################################################################################
 Input:
    frames      - dcd file
    pdb_list    - list with pdb file lines
    psf_list    - list with pdb file lines
    upper_types - type of top atoms, in this case Nitrogen
    mid_type    - type of middle type atom that connects paba with the cellulose
 ###################################################################################
 Description:
  Iterate over each trajectory and over each line in the pdb to find the number ids
  of the selected atoms. This ids are based in the simulation, not on the list indexes.
  In this case the ids are 4064 and 4076.
  uppTypePos = 4064 = ['4065', 'MAIN', 'A', 'PAB', 'N', 'N000', '-0.551000', '14.0067', '0']
  midTypePos = 4076 = ['4077', 'MAIN', 'A', 'PAB', 'C', 'C150', '+0.171900', '12.0107', '0']
 ###################################################################################
  Output:
    midTypePos - list with the number ids of middle type atoms
    uppTypePos - list with the number ids of upper type atoms
    contPABA   - number of pabas in one frame
 ###################################################################################
"""
def getIndexUppAndMidTypes(frames, pdb_list, psf_list, upper_types, mid_type):
    contPABA = 0
    midTypePos = list()  # list to store indexs of mid_types atoms
    uppTypePos = list()  # list to store indexs of upper_types atoms
    for frame in frames: # iterate through dcd file
        for atom in range(len(pdb_list)): # for range in the amount of atoms
            try:
                # if atom is PABA go for it
                if psf_list[atom][3] == 'PAB':
                    # if atom is one of the Nitrogen or a base type atom go for it and store it
                    if (psf_list[atom][5])[:len(upper_types)] in upper_types:
                        uppTypePos.append(int(pdb_list[atom][1])) # store atom number id from simulation
                    if (psf_list[atom][5])[:len(mid_type)] == mid_type:
                        contPABA += 1
                        midTypePos.append(int(pdb_list[atom][1])) # store atom number id from simulation
            except IndexError:
                pass
        break
    return midTypePos, uppTypePos, contPABA
#------------------------------------------------------------------------
"""
 ###################################################################################
 Function: getIndexUppAndMidTypes(frames, pdb_list, psf_list, upper_types, mid_type)
 ###################################################################################
 Input:
    frames      - dcd file
    pdb_list    - list with pdb file lines
    psf_list    - list with pdb file lines
    upper_types - type of top atoms, in this case Nitrogen
    mid_type    - type of middle type atom that connects paba with the cellulose
 ###################################################################################
 Description:
  Iterate over each trajectory and over each line in the pdb to find the number ids
  of the selected atoms. This ids are based in the simulation, not on the list indexes.
  In this case the ids are 4064 and 4076.
  uppTypePos = 4064 = ['4065', 'MAIN', 'A', 'PAB', 'N', 'N000', '-0.551000', '14.0067', '0']
  midTypePos = 4076 = ['4077', 'MAIN', 'A', 'PAB', 'C', 'C150', '+0.171900', '12.0107', '0']
 ###################################################################################
  Output:
    midTypePos - list with the number ids of middle type atoms
    uppTypePos - list with the number ids of upper type atoms
    contPABA   - number of pabas in one frame
 ###################################################################################
"""
def getUpperHeights(frames, uppPos):
    heights   = []
    z_coor    = []
    contUpper = 0
    for frame in frames:
        for i in uppPos:
            #print frame[0][i-1],frame[1][i-1],frame[2][i-1]
            z_coor.append(float(frame[2][i-1])) # append to list
            contUpper += 1 # add to the counter
            if contUpper == 3:
                heights.append(np.array(heapq.nlargest(2,z_coor)).mean()) # print mean to store it in a file
                z_coor    = [] # empty list
                contUpper = 0  # renew counter
        #print "End frame"
    return heights
#------------------------------------------------------------------------
"""
 ###################################################################################
 Function: getIndexUppAndMidTypes(frames, pdb_list, psf_list, upper_types, mid_type)
 ###################################################################################
 Input:
    frames      - dcd file
    pdb_list    - list with pdb file lines
    psf_list    - list with pdb file lines
    upper_types - type of top atoms, in this case Nitrogen
    mid_type    - type of middle type atom that connects paba with the cellulose
 ###################################################################################
 Description:
  Iterate over each trajectory and over each line in the pdb to find the number ids
  of the selected atoms. This ids are based in the simulation, not on the list indexes.
  In this case the ids are 4064 and 4076.
  uppTypePos = 4064 = ['4065', 'MAIN', 'A', 'PAB', 'N', 'N000', '-0.551000', '14.0067', '0']
  midTypePos = 4076 = ['4077', 'MAIN', 'A', 'PAB', 'C', 'C150', '+0.171900', '12.0107', '0']
 ###################################################################################
  Output:
    midTypePos - list with the number ids of middle type atoms
    uppTypePos - list with the number ids of upper type atoms
    contPABA   - number of pabas in one frame
 ###################################################################################
"""
def getMidTypePos(frames, midTypePos, pdb_list, psf_list, base_type):
    midTypeList = list()
    for frame in frames:
        for i in midTypePos:
            temp    = 100 # random max value - cambiar por el maximo de python
            tempPos = 0
            for atom in range(len(pdb_list)):
                try:
                    if psf_list[atom][3] == 'CEL' and (psf_list[atom][5])[:len(base_type)] == base_type:
                        dist = math.sqrt((frame[0][i-1] - frame[0][atom])**2 + (frame[1][i-1] - frame[1][atom])**2 + (frame[2][i-1] - frame[2][atom])**2)
                        if dist < temp:
                            temp    = dist
                            tempPos = int(pdb_list[atom][1])
                except IndexError:
                    pass
            midTypeList.append(tempPos)
        break # just run for one frame, positions will be the same for each frame
    return midTypeList
#------------------------------------------------------------------------
"""
 ###################################################################################
 Function: getIndexUppAndMidTypes(frames, pdb_list, psf_list, upper_types, mid_type)
 ###################################################################################
 Input:
    frames      - dcd file
    pdb_list    - list with pdb file lines
    psf_list    - list with pdb file lines
    upper_types - type of top atoms, in this case Nitrogen
    mid_type    - type of middle type atom that connects paba with the cellulose
 ###################################################################################
 Description:
  Iterate over each trajectory and over each line in the pdb to find the number ids
  of the selected atoms. This ids are based in the simulation, not on the list indexes.
  In this case the ids are 4064 and 4076.
  uppTypePos = 4064 = ['4065', 'MAIN', 'A', 'PAB', 'N', 'N000', '-0.551000', '14.0067', '0']
  midTypePos = 4076 = ['4077', 'MAIN', 'A', 'PAB', 'C', 'C150', '+0.171900', '12.0107', '0']
 ###################################################################################
  Output:
    midTypePos - list with the number ids of middle type atoms
    uppTypePos - list with the number ids of upper type atoms
    contPABA   - number of pabas in one frame
 ###################################################################################
"""
def getBaseCoordinates(frames, midTypeList):
    baseTypeList = list()
    contTot  = 0
    #print midTypeList
    for frame in frames:
        for atom in midTypeList:
            contTot += float(frame[2][atom-1])
            baseTypeList.append(float(frame[2][atom-1]))
    return baseTypeList, contTot
#------------------------------------------------------------------------
"""
 ###################################################################################
 Function: calculateAverageHeight(total, contPABA, heights)
 ###################################################################################
 Input:
   total    - sum of all the heights together
   contPABA - number of pabas in a frame
   heights  - heights of each paba through all the frames
 ###################################################################################
 Description:
   Calculate the average of all the heights and substract that to the final height. 
   Results are outputted to the terminal.
 ###################################################################################
 Function: calculateNormalizedHeights(heigths, baseHeights)
 ###################################################################################
  Input:
   heights     - heights of each paba through all the frames
   baseHeights - heights of each base atom through all the frames
 ###################################################################################
 Description:
   Calculate the height minus the base atom to get a normalized height of the ligand.
   Results are outputted to the terminal.
 ###################################################################################
"""
def calculateAverageHeight(total, contPABA, heights):
    baseHeight = float(total) / (float(contPABA) * 1000.0)
    for h in heights:
        print h - baseHeight
#------------------------------------------------------------------------
def calculateNormalizedHeights(heigths, baseHeights):
    for h in range(len(heigths)):
        print heigths[h] - baseHeights[h] 
#------------------------------------------------------------------------   





midTypePos, uppTypePos, pabaNum = getIndexUppAndMidTypes(dcd[0], pdb_list, psf_list, upper_types, mid_type)
upperHeights                    = getUpperHeights       (dcd[1], uppTypePos)
midTypeList                     = getMidTypePos         (dcd[2], midTypePos, pdb_list, psf_list, base_type)
baseHeights, baseTot            = getBaseCoordinates    (dcd[3], midTypeList)
calculateNormalizedHeights(upperHeights, baseHeights)
#calculateAverageHeight(baseTot, pabaNum, upperHeights)

