#!/bin/bash
#-------------------------------------------------------------------------------------------------------
# Script developed to extract certain frames from the dcd file taken from a simulation. This script
# uses catdcd to cut a frame from the dcd and convert it into a pdb file. This pdb file has the
# coordinates of the atoms and lets the user calculate variations like height, movement, and frequency.
# For more information on catdcd please visit http://www.ks.uiuc.edu/Development/MDTools/catdcd/
# This program assumes that catdcd will be installed in your system. If it is not installed, you can
# change the path in the command to where it is located.
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
#-------------------------------------------------------------------------------------------------------

# Enter the names of the dcd and the initial pdb file
# catdcd - path to catdcd binary file E.g.: /Desktop/Cellulose/SimulationsDone/catdcd
# if it is installed it can be plain catdcd
# path - path for the resulting pdb files E.g.: /Desktop/Cellulose/SimulationsDone/Last1000PDBs/
# dcd_file - path to dcd file E.g.: /Desktop/Cellulose/SimulationsDone/name.dcd
# pdb_file - path to pdb file E.g.: /Desktop/Cellulose/SimulationsDone/name.pdb
catdcd="None"
path="None"
dcd_file="None"
pdb_file="None"

# Simulation name (descriptive name that will be in the output files header)
# E.g. Pentane100
description="None"

# Select the initial frame from the dcd file
InitFrame=9000
LastFrame=1000

# Enter in a loop to save the pdb file from each frame selected
for (( i=1; i<=$LastFrame; i++ ))
do
    # Update the frame with a counter
    frame=$[$i+$InitFrame]

    # Run the catdcd tool with the given parameters
    $catdcd -o $path${description}_Frame${frame}.pdb -otype pdb -s $pdb_file -stype pdb -first $frame -last $frame $dcd_file

    # This will produce a bunch of pdb files in the current path representing the selected frame

done
