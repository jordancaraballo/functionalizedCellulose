"""
Purpose: Take spacer arm coordinates not including hydrogens or the para-aminobenzamidine ring.
Script is written in python2.7 for compatibility with Wolffia.
To run script script:
    python2 ./spacerArmCoordinates.py mixture_name working_path
    python2 ./spacerArmCoordinates.py Tetradecane005 /home/fulano/Tetradecane
    Note that the mixture is going to be automatically appended at the end of the path.
"""
#!/bin/python
#-------------------------------------------------------------------------------------------------------
# Imports
import sys, os                            # environmental variables
import heapq                              # import for median
import math                               # import for distance calculation
import numpy as np                        # include numpy
import csv
#-------------------------------------------------------------------------------------------------------
# Additional Imports
sys.path.append('/home/jordancaraballo/WolffiaApp/wolffia/')  # add wolffia repository for DCDReader utility
from lib.io.dcd import DCDReader                              # include wolffia DCDReader library
#-------------------------------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------------------------------
# Reading and parsing files
def readPDB(pdbFile):
    return [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:]
def readPSF(psfFile):
    return [line.strip().split() for line in open(psfFile,'r').readlines() ][6:]
def readDCD(dcdFile):
    return DCDReader(dcdFile)
# Cleaning and parsing coordinates from frames
def getCoordinatesSingleFrameSpacerArm(pdb_list, psf_list, singleDCD, typesSA):
    atoms_xyz_dict = dict()
    molecules_list = list()
    for frame in singleDCD:
        for atom in range(len(pdb_list)): # iterate over every atom in the frame
            try:
                if psf_list[atom][3] == 'PAB' and (psf_list[atom][5])[:3] in typesSA: # if atom is from spacer arms
                    #print "Atom type: ", (psf_list[atom][5])[:3], "XCoor: ", float(frame[0][int(pdb_list[atom][1])-1])
                    #print "Atom type: ", (psf_list[atom][5])[:3], "XCoor: ", float(frame[0][int(pdb_list[atom][1])-1])

                    atoms_xyz_dict[(psf_list[atom][5])[:3]] = (float(frame[0][int(pdb_list[atom][1])-1]),
                                                               float(frame[1][int(pdb_list[atom][1])-1]),
                                                               float(frame[2][int(pdb_list[atom][1])-1]))
                if len(atoms_xyz_dict) == len(typesSA):
                    print len(atoms_xyz_dict)
                    molecules_list.append(atoms_xyz_dict)
                    atoms_xyz_dict.clear()

            except IndexError: # some lines are trash lines from the pdb file
                pass
    return molecules_list

# Get coordinates from a single frame, just the spacer arm atoms. No hydrogens, no nothing.
def getCoordinatesSingleFrameSpacerArmOld(pdb_list, psf_list, singleDCD, molLength):
    tempCoors   = list() # list to store temporary coordinates for each pABA
    coordinates = list()
    tempcoordinates_types = list() # list to store final coordinates types
    coordinates_types     = list()
    types_spacerarm = ['C60', 'C70','C80','C90','C10','C11', 'C12','O20','C13','C14','C15'] # types of spacer arm

    for frame in singleDCD: # iterate over every frame
        for atom in range(len(pdb_list)): # iterate over every atom in the frame
            try:
                if psf_list[atom][3] == 'PAB' and (psf_list[atom][5])[:3] in types_spacerarm: # if atom is from spacer arms
                    tempCoors.append(float(frame[0][int(pdb_list[atom][1])-1])) # append x coordinates
                    tempCoors.append(float(frame[1][int(pdb_list[atom][1])-1])) # append y coordinates
                    tempCoors.append(float(frame[2][int(pdb_list[atom][1])-1])) # append z coordinates
                    tempcoordinates_types.append(str((psf_list[atom][5])[:3]))

                    if len(tempCoors) == len(types_spacerarm)*3: # if tempCoors has all the atoms from the molecule
                        coordinates.append(tempCoors) # send atoms to the main list
                        coordinates_types.append(tempcoordinates_types)
                        tempCoors = []                # empty temp list
                        tempcoordinates_types = []
            except IndexError: # some lines are trash lines from the pdb file
                pass
    return coordinates, coordinates_types # multidimensional array

#---------------------------------------------------------------------------
# Main
#---------------------------------------------------------------------------
if __name__ == "__main__":
    # Global variables
    mix = sys.argv[1]                            # mixture name
    workPath = sys.argv[2] + mix + "/"           # path to mixture - parental directory
    psfFile  = workPath + mix + ".psf"           # name of the psf file
    pdbFile  = workPath + mix + "_lastFrame.pdb" # header of the pdb files
    dcdFile  = workPath + mix + "_LastFrame.dcd" # dcd file

    # Reading and cleaning files
    pdb_list = readPDB(pdbFile)
    psf_list = readPSF(psfFile)
    dcd_list = readDCD(dcdFile)

    # Processing atoms
    types_spacerarm = ['C60', 'C70','C80','C90','C10','C11', 'C12','O20','C13','C14','C15'] # types of spacer arm
    atoms = getCoordinatesSingleFrameSpacerArm(pdb_list, psf_list, dcd_list, types_spacerarm)
    print len(atoms)
