#!/bin/python
#-------------------------------------------------------------------------------------------------------
# Script developed to functionalize cellulose nanocrystals with pABA ligands working as spacer arms.
# This script uses functionalities from the wolffia project that can be found at
# https://github.com/compMathUPRH/wolffia. In order to run this script you will need to clone wolffia's
# project first and then add it to your path as a python path variable. The aim of this project is to study
# the frequencies of ligands over this cellulose surfaces at different densities.
# The output of this program is a set of .wfm files that will have the configurations for each
# simulation for NAMD. You can easily see the mixture in wolffia. This script is directed to 5-atoms spacer
# arms made of para-aminobenzamidine.
# Authors: Jordan Alexis Caraballo-Vega  University of Puerto Rico at Humacao
#          Jose O. Sotero-Esteva         University of Puerto Rico at Humacao
#-------------------------------------------------------------------------------------------------------

# Import libraries
from lib.chemicalGraph.molecule.polymer.Cellulose import Cellulose
from lib.chemicalGraph.molecule.solvent.PABA import PABA_Tetradecane
from lib.chemicalGraph.Mixture import Mixture
import random
random.seed()
from math import ceil

# Elements for crystal construction - this options are the ones that initialize the crystal arrangement
POLY_LENGTH  = 4 # number of the polymer chains in z axis (from lower to top)
HORIZONTAL_N = 6 # number of the polymer chains in x axis (from left to right)
VERTICAL_N   = 4 # number of chains in y axis (from front to the back)

# Elements for crystal arrangement
dx  = 8.0 # number of armstrongs to move crystals in x axis
dy  = 0.0 # number of armstrongs to move crystals in y axis
dz  = 4.  # number of armstrongs to move crystals in z axis
dxz = 3   # common movement to use in a loop

# Lists of densities and and names for output files
densities = [100.0,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05]
d_names   = ["100","09","08","07","06","05","04","03","02","01","005"]

# Iterate through every density
for i in range(len(densities)):

	# Density value
	PABAratio = densities[i]
	FileName  = "Cellulose_PABA_Tetradecane_D" + d_names[i] + ".wfm"

	# Define mixture
	mixture = Mixture()

	# Generate crystal minus last layer
	monomer = Cellulose(POLY_LENGTH)
	for j in range(VERTICAL_N-1):
	  for i in range(HORIZONTAL_N):
	    newMonomer = monomer.copy()
	    newMonomer.moveBy([dx * i + dxz * (j%2), dy * (i%2), dz * j])
	    mixture.add(newMonomer)
	    print "layer " + str(j) + ", row " + str(i)

	# Generate last layer
	for i in range(HORIZONTAL_N):
	  newMonomer = monomer.copy()
	  newMonomer.moveBy([dx * i + dxz * ((VERTICAL_N-1)%2), dy * (i%2), dz * (VERTICAL_N-1)])
	  mixture.add(newMonomer)
	  #print "last layer, row " + str(i)

	#find OHs at the Cellulose nanocrystal
	OHs = list()
	for m in mixture.moleculeGenerator():
	  for atom in m:
	    t = m.getAtomAttributes(atom).getInfo().getType()
	    if t == "H13" or t == "H19":
	     oxygen = m.neighbors(atom)[0]
	     coords = m.getAtomAttributes(oxygen).getCoord()
	     if coords[2] > (VERTICAL_N-1) * dz:  # get only the top O's
	       OHs.append( (m,atom) )

	# Print number of detected OHs
	print "Total Active OHs ", len(OHs)

	#adding PABA
	paba = PABA_Tetradecane()
	paba.rename("PABAtetradecane")
	paba.removeAtom(327)  #remove hydrogen
	paba.removeAtom(241)  #remove oxygen
	paba.moveBy([-x for x in paba.center()])
	paba.rotateDeg(90.0, 0.0, 0.0)
	paba.moveBy([0.0,0.0,12.0])
	chrgToSubstract = 0.4181  # computed with Wolffia
	#chrgToSubstract = 0.208801
	chrgToSubstractO200 = 0.0181  # computed with Wolffia

	# select functionalization sites
	ohPlaces = range(len(OHs))
	for i in range(int(ceil(len(ohPlaces)*(1-PABAratio)))):
	  del ohPlaces[random.randint(0,len(ohPlaces))-1]
	print len(ohPlaces), "OHs selected for functionalization"

	for i in ohPlaces:
	  oh = OHs[i]
	  oxygen = oh[0].neighbors(oh[1])[0]
	  oh[0].removeAtom(oh[1])  # remove hydrogen from celulose
	  #print "PABA Ox " + str(oxygen)
	  coords = oh[0].getAtomAttributes(oxygen).getCoord()
	  atype = oh[0].getAtomAttributes(oxygen).getInfo().getType()
	  if coords[2] > (VERTICAL_N-1) * dz:  # get only the top O's
	    newPABA = paba.copy()
	    newPABA.rotateDeg(0.0, 0.0, 180.0 * random.random())
	    newPABA.moveBy(coords)

	    #state.addMolecule(newPABA)
	    oh[0].merge(newPABA, [(oxygen,114)])

	    #neutralize total charge
	    ch = oh[0].atom_attributes(oxygen).getInfo().getCharge()
	    if atype == "O2":
	      oh[0].atom_attributes(oxygen).getInfo().setCharge(ch-(-0.208801))
	    else:
	      oh[0].atom_attributes(oxygen).getInfo().setCharge(ch-(-0.208801))

        print "Force field ", oh[0].getForceField()._BONDS

        oh[0].getForceField().setBond(('C15', 'O6'), 1.42, 1)
        oh[0].getForceField().setBond(('C15', 'O6'), 428, 0)
        oh[0].getForceField().setBond(('C15', 'O2'), 1.42, 1)
        oh[0].getForceField().setBond(('C15', 'O2'), 428, 0)

        oh[0].getForceField().setAngle(('C15', 'O2', 'C2'), 131.929, 0)
        oh[0].getForceField().setAngle(('C15', 'O2', 'C2'), 78., 1)
        oh[0].getForceField().setAngle(('C15', 'O6', 'C6'), 131.929, 0)
        oh[0].getForceField().setAngle(('C15', 'O6', 'C6'), 104., 1)

        oh[0].getForceField().setAngle(('C14', 'C15', 'O6'), 113.391, 0)
        oh[0].getForceField().setAngle(('C14', 'C15', 'O6'), 40., 1)
        oh[0].getForceField().setAngle(('C14', 'C15', 'O2'), 113.391, 0)
        oh[0].getForceField().setAngle(('C14', 'C15', 'O2'), 40., 1)

	mixture.save(FileName)
