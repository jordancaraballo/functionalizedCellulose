import networkx as nx
import matplotlib.pyplot as plt

# Script to parse and modify PAB. Designed for python2 since wolffia is python2.
import sys, os                                     # environmental variables
import itertools                                   # iterate over list
import heapq                                       # import for median
import math                                        # import for distance calculation
import numpy as np                                 # include numpy
import subprocess                                  # execute bash commands if needed
#----------------------------------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------------------------------
def parsePSF_toList(psfFile): # input: string with psf file name, output: list of file lines splitted
    print "Parsing PSF..."
    return [line.strip().split() for line in open(psfFile,'r').readlines() ][6:] # list with psf file lines
#----------------------------------------------------------------------------------------------------------------
def parsePDB_toList(pdbFile): # input: string with pdb file name, output: list of file lines splitted
    print "Parsing PDB..."
    return [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:] # list with pdb file lines
#----------------------------------------------------------------------------------------------------------------
def parseDCD_toList(dcdFile): # input: string with dcd file name, output: list of dcds copies
    print "Parsing DCD..."
    return [DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile)] # list with dcd files
#----------------------------------------------------------------------------------------------------------------
# combination of pdb and psf file to parse paba coordinates
# input: dcd, pdb list, psf list, residue to match. output: list of atoms
def getAtoms(psf_list, residue, tcolumn=5, mcolumn=3, mol="PAB"):
    print "Extracting Atoms for residue " + residue + "In mol " + mol + "..."
    atomsInfo = list() # list to store atoms
    for atom in psf_list: # for range in the number of atoms
        try:
            if atom[tcolumn][:len(residue)] == residue and atom[mcolumn][:len(mol)] == mol: # if atom is PABA go for it
                atomsInfo.append(int(atom[0])) # append atom line to list
        except IndexError:
            pass
    print "Total Atoms: ", len(atomsInfo)
    return atomsInfo
#----------------------------------------------------------------------------------------------------------------
# combination of pdb and psf file to parse paba coordinates
# input: dcd, pdb list, psf list, residue to match. output: list of atoms
def getAtomsMol(psf_list, mol, mcolumn=3):
    print "Extracting Atoms for molecule " + mol + "..."
    atomsInfo = list() # list to store atoms
    for atom in psf_list: # for range in the number of atoms
        try:
            if atom[mcolumn][:len(mol)] == mol: # if atom is PABA go for it
                atomsInfo.append(int(atom[0])) # append atom line to list
        except IndexError:
            pass
    print "Total Atoms: ", len(atomsInfo)
    return atomsInfo
#----------------------------------------------------------------------------------------------------------------
# Note: to find numBonds you can type cat /path/Tetradecane05.psf | grep NBOND.
def getBonds(psf_list, numBonds): # input: psf list, number of bonds in file. output: list with bond connections
    lst = psf_list[psf_list.index([numBonds, '!NBOND'])+1:]
    lst = lst[:lst.index([])]
    pairs = []
    for line in lst:
        for i in range(0,len(line), 2):
            pairs.append((int(line[i]), int(line[i+1])))
    return pairs
#----------------------------------------------------------------------------------------------------------------
def getSpacerArmGraph(graph, baseAtoms, upperAtoms):
    sparmList = list()
    for paba in graph:
        for start in  baseAtoms:
            for end in upperAtoms:
                try:
                    sparm = nx.shortest_path(paba, source=start, target=end)
                    sparmList.append(sparm)
                except: pass
    return sparmList
#----------------------------------------------------------------------------------------------------------------


# Main
if __name__ == "__main__":
    
    psfFile = "/home/jordancaraballo/Documents/Research/Cellulose/Tetradecane/Tetradecane100/Tetradecane100.psf" # name of the psf file
    psfList = parsePSF_toList(psfFile) # takes the list of psf lines

    pabaAtoms = getAtomsMol(psfList, 'PAB') # takes paba atoms from psflist
    
    # Tetradecane - Base: C15, Upper: C60
    baseAtoms  = getAtoms(psfList, 'C15', tcolumn=5) # takes C15 (base carbon) atom from each molecule
    upperAtoms = getAtoms(psfList, 'C60', tcolumn=5) # takes C06 (upper carbon) atom from each molecule

    numBonds  = (subprocess.check_output("cat " + psfFile + "| grep NBOND | cut -d \" \" -f4", shell=True))[:-1]
    bondConnects = getBonds(psfList, numBonds) # get bond connections list
    
    G = nx.Graph()
    G.add_edges_from(bondConnects)

    Gsub = G.subgraph(pabaAtoms)
    #nx.draw(Gsub)
    
    comp = list(Gsub.subgraph(c) for c in nx.connected_components(Gsub))
    print getSpacerArmGraph(comp,baseAtoms,upperAtoms)
#for paba in comp:
#    for comienzo in  c15Atoms:
#       for final in c06Atoms:
#           try:
#               pelo = nx.shortest_path(paba, source=comienzo, target=final)
#              print(pelo)
#           except: pelo = "(nada)"
 
#nx.draw(comp[0])
