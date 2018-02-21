#!/bin/python
"""
#-------------------------------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------------------------
"""
import heapq
import math
import sys, os
sys.path.append('/home/Fulano/wolffia/')
from lib.io.dcd import DCDReader
import numpy as np

# Global variables
mix = "Pentane100" # number of the mixture - e.g Pentane09
workPath = "/home/Fulano/" + mix + "/"
psfFile  = workPath + mix + ".psf"                 # name of the psf file
pdbFile  = workPath + mix + "_lastFrame.pdb"       # header of the pdb files
dcdFile  = workPath + mix + "_Last1000Frames.dcd"  # dcd file

# Atom types to detect
# Pentane mid_type = C8, Tetradecane mid_type = C15
upper_types = ["N0", "N1"] # atoms to get coordinates that we want
base_type   = 'O6'         # atom to calculate the height of pABA
mid_type    = 'C8'         # atom that connects pABA and Cellulose

# Creates list of each line on the psf file - this includes each atom info
psf_list = [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]
pdb_list = [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:]
dcd  = [DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile)]
#-------------------------------------------------------------------------------------
# Search for PABA two top atoms (two external Nitrogens) and get a mean of them.
# Then do a dendogram of each of them that would express their height.
# The coor to get is in the z axis.
#--------------------------------------------------------------------------------------
z_coor      = list() # temporary array to calculate coordinates
heights     = list() # list to store heights from pABA's
midTypeList = list() # list to store indexs of mid_types atoms
baseHeights = list()
contPABA    = 0

#------------------------------------------------------------------------
midTypePos = list() # list to store indexs of mid_types atoms
uppTypePos = list() # list to store indexs of upper_types atoms
for frame in dcd[0]:
    for atom in range(len(pdb_list)):
        try:
            if psf_list[atom][3] == 'PAB':
                if (psf_list[atom][5])[:2] in upper_types:
                    uppTypePos.append(atom)
                if (psf_list[atom][5])[:2] == mid_type:
                    contPABA += 1
                    midTypePos.append(atom)
        except IndexError:
            pass
    break
#------------------------------------------------------------------------
z_coor    = []
contUpper = 0
for frame in dcd[1]:
    for i in uppTypePos:
        z_coor.append(float(frame[2][i])) # append to list
        contUpper += 1 # add to the counter
        if contUpper == 3:
            heights.append(np.array(heapq.nlargest(2,z_coor)).mean()) # print mean to store it in a file
            z_coor    = [] # empty list
            contUpper = 0  # renew counter
#------------------------------------------------------------------------
midTypeList = list()
for frame in dcd[2]:
    for i in midTypePos:
        temp    = 100
        tempPos = 0
        for atom in range(len(pdb_list)):
            try:
                if psf_list[atom][3] == 'CEL' and (psf_list[atom][5])[:2] == base_type:
                    dist = math.sqrt((frame[0][i] - frame[0][atom])**2 + (frame[1][i] - frame[1][atom])**2 + (frame[2][i] - frame[2][atom])**2)
                    if dist < temp:
                        temp    = dist
                        tempPos = atom
            except IndexError:
                pass

        midTypeList.append(tempPos)
    break
#------------------------------------------------------------------------
baseTypeList = list()
contPABA = 0
contTot  = 0
for frame in dcd[3]:
    for atom in midTypeList:
        contPABA += 1
        contTot += float(frame[2][atom])
        baseTypeList.append(float(frame[2][atom]))

baseHeight = contTot / float(contPABA)
for h in heights:
    print h - baseHeight
