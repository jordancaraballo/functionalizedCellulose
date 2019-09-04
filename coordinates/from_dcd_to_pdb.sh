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

## Path to catdcd file - default is assuming it is installed
catdcd="catdcd"

## Iterate to get arguments from command line
## Ignoring input validation for now since this is called from another file
while [ -n "$1" ]; do # while loop starts

    case "$1" in

    --dcd) # input dcd file with simulation frames
        dcd="$2"
        shift
        ;;
    --pdb) # input pdb file with simulation atoms, this will facilitate the conversion of pdb frames
        pdb="$2"
        shift
        ;;
    --first) # first frame to start extracting coordinates
        first="$2"
        shift
        ;;
    --last) # last frame to start extracting coordinates
        last="$2"
        shift
        ;;
    --out) # output description
        out="$2"
        shift

        break
        ;;

    *) echo "Option $1 not recognized" ;;

    esac
    shift

done

## Enter in a loop to save the pdb file from each frame selected
for (( i=1; i<=$last; i++ ))
do
    frame=$[$i+$first] # Update the frame with a counter
    $catdcd -o ${out}_Frame${frame}.pdb -otype pdb -s $pdb -stype pdb -first $frame -last $frame $dcd
done