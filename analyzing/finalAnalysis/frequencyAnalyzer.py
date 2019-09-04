# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 08:09:47 2019

@author: jordancaraballo
"""
# Script to parse and modify PAB. Designed for python2 since wolffia is python2.
# Authors: Jordan A. Caraballo Vega, Jose O. Sotero Esteva
#----------------------------------------------------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------------------------------------------------
import numpy as np                           
from AnalyzerFunctions import parsePSF_toList, parsePDB_toList
from AnalyzerFunctions import getAtoms
from AnalyzerFunctions import getNitrogenCoordinates
from AnalyzerFunctions import frequencyChart
#----------------------------------------------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    system       = "Pentane" # "Heptane" "TetradecaneBig" "Tetradecane"
    base_atom    = 'C80' # C15 Tetradecane, C10 Heptane, C80 Pentane
    densities    = ["005","01","02","03","04","05","06","07","08","09","100"] 
    base_colors  = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','cyan', 'fuchsia', "orange"]
    #base_colors = ['blue','green','red','cyan','magenta','yellow','black','indigo','burlywood','darksalmon','darkviolet']
    
    finalSet     = list() # final list with the set of coordinates
    colors       = list() # final list with the set of colors per concentration
    minList      = list()
    maxList      = list()

    for den in densities:
        
        print "---------"
        ### Getting files information and loading data into lists
        psfFile = "data/" + system + "/" + system + den + ".psf"
        psfList = parsePSF_toList(psfFile) # takes the list of psf lines
        
        pdbFile = "data/" + system + "/" + system + den + "_LastFrame.pdb"
        pdbList = parsePDB_toList(pdbFile)

        ### Starting to define spacer arms information and coordinates
        baseAtoms = getAtoms(psfList, base_atom, tcolumn=5) # takes C15 (base carbon) atom from each molecule
        N1Atoms   = getAtoms(psfList, 'N10', tcolumn=5) # takes N1 (second nitrogen) atom from each molecule
        
        ### Substract heigth minus base
        atomCoordinatesN1  = getNitrogenCoordinates(N1Atoms, psfList, pdbList) # get atoms coordinates (np array)
        atomCoordinatesC15 = getNitrogenCoordinates(baseAtoms, psfList, pdbList) # get atoms coordinates (np array)
        atomCoordinates    = np.array(atomCoordinatesN1) - np.array(atomCoordinatesC15)
        
        atomCoordinates    = atomCoordinatesN1
        
        minList.append(min(atomCoordinates))
        maxList.append(max(atomCoordinates))
        finalSet.append(atomCoordinates)
        
    finalSet = np.array(finalSet)   
    dMin = min(minList)
    dMax = max(maxList)
    
    frequencyChart(finalSet, dMin, dMax, base_colors, system)
     

"""       
        ### Starting to define graph and shortests path
        numBonds  = (subprocess.check_output("cat " + psfFile + "| grep NBOND | cut -d \" \" -f4", shell=True))[:-1]
        bondConnects = getBonds(psfList, numBonds) # get bond connections list
    
        G = nx.Graph()
        G.add_edges_from(bondConnects)
        Gsub = G.subgraph(pabaAtoms)
        comp = list(Gsub.subgraph(c) for c in nx.connected_components(Gsub))
        
        ### Starting to get atom indices and coordinates from the spacer arm only
        graphIndexes = getSpacerArmGraph(comp,baseAtoms,upperAtoms) # list of atom indexes in psf file (list)
        atomCoordinates = getPABACoordinates(graphIndexes, psfList, pdbList) # get atoms coordinates (np array)

        colors.extend([base_colors[densities.index(den)]] * atomCoordinates.shape[0]) # create list of colors        
        rotatedCoordinates = rotate_matrix_full(atomCoordinates) # Rotate coordinates with respect to x, y, and z axis
        reshapedCoordinates = reshapeCoordinates(rotatedCoordinates) # Reshape coordinates for additional visualizations
        finalSet.append(reshapedCoordinates) # Append new atoms to greater list
        
    ### Append coordinates to final set of atoms and visualize
    finalSet = np.array(finalSet) # convert to numpy array (from multi dim. list)
    finalSet = np.concatenate(np.array(finalSet),axis=0) # concatenate the multi dim into a set of one array
    colors = np.array(colors) # convert colors list into numpy array
    
    ### Classification algorithms - See shapeAnalyzerFunctions.py for additional classifiers
    # kmeans(finalSet, colors), pca2d(finalSet, colors), pca2d_OnClick(finalSet, colors)
    # pca3d(finalSet, colors), hierarchical(finalSet), tSNE(finalSet, colors, 500)
    pca2d_OnClickFit(finalSet, colors)
"""
