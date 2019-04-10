# Script to parse and modify PAB. Designed for python2 since wolffia is python2.
# Authors: Jordan A. Caraballo Vega, Jose O. Sotero Esteva
#----------------------------------------------------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------------------------------------------------
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys, os                                     # environmental variables
import itertools                                   # iterate over list
import heapq                                       # import for median
import math                                        # import for distance calculation
import numpy as np                                 # include numpy
import subprocess                                  # execute bash commands if needed
from scipy import optimize                         # for least square 

from sklearn import datasets
from sklearn.cluster import KMeans
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
    print "Extracting Atoms for residue " + residue + " In mol " + mol + "..."
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
#
def getPABACoordinates(indexesList, psfList, pdbList):
    atomsCoordinates = list()
    for paba in indexesList:
        temp = list()
        for atom in paba:
            temp.append([float(pdbList[atom-1][6]), float(pdbList[atom-1][7]), float(pdbList[atom-1][8])])
            #print "Atom Index: ", atom, " PDB Atom: ", pdbList[atom-1][6], pdbList[atom-1][7], pdbList[atom-1][8]
        atomsCoordinates.append(temp)
    return np.array(atomsCoordinates)
        
#----------------------------------------------------------------------------------------------------------------
# Purpose: translate coordinates to the origin based on lowest pABA atom
def translate_coordinates_lowestAtom(molecule_array, index_base=0):
    return molecule_array - molecule_array[index_base]
#---------------------------------------------------------------------------
# Purpose: rotate coordinates over RyRx axis
def rotate_matrix_RyRxHighest(molecule, index_highest=-1):
    V1,V2,V3 = molecule[index_highest] # coordinates defined as vector representations

    # Solving for rotational matrix: a = cos x, b = sin x, c = cos x, d = sin x
    a = math.sqrt(V3**2 / (V2**2 + V3**2))
    b = math.sqrt(1 - V3**2 / (V2**2 + V3**2))
    c = math.sqrt(1 - V1**2 / ((V3 * a + V2 * b)**2 + V1**2))
    d = math.sqrt(V1**2 / ((V3 * a + V2 * b)**2 + V1**2))

    return molecule.dot(np.array([[c, d * b, d * a], [0, a, -b], [-d, b * c, a * c]])) # return rotated molecule matrix
#---------------------------------------------------------------------------
# Purpose: rotate coordinates over Rz axis
def rotate_matrix_RzMiddle(molecule, index_middle=5):
    V1,V2,V3 = molecule[index_middle] # coordinates defined as vector representations

    # Solving for rotational matrix: e = cos z, f = sin z
    e = math.sqrt(V2**2 / (V1**2 + V2**2))
    f = math.sqrt(1 - V2**2 / (V1**2 + V2**2))

    return molecule.dot(np.array([[e, -f, 0], [f, e, 0], [0, 0, 1]])) # return rotated molecule matrix
#---------------------------------------------------------------------------
# Purpose: Find A of the spiral equation = A: height of the upper atom
def getVarSpiral_A(molecule, index_highest=-1, index_base=0):
    return molecule[index_highest][2] - molecule[index_base][2]
#---------------------------------------------------------------------------
# Playing with KMeans Algorithm
def kmeans(atoms):

    # Loading dataset
    #iris_df = datasets.load_iris()
    #print iris_df

    # Declaring Model
    model = KMeans(n_clusters=3)
    
    # Fitting Model
    #model.fit(iris_df.data)
    model.fit(atoms)
    
    # Predicitng a single input
    #predicted_label = model.predict([[7.2, 3.5, 0.8, 1.6]])
    
    # Prediction on the entire data
    #all_predictions = model.predict(iris_df.data)
    all_predictions = model.predict(atoms)
    
    # Printing Predictions
    #print(predicted_label)
    print(all_predictions)

def hierarchical(atoms):
    # Importing Modules
    from scipy.cluster.hierarchy import linkage, dendrogram
    import matplotlib.pyplot as plt

    """
    Perform hierarchical clustering on samples using the
    linkage() function with the method='complete' keyword argument.
    Assign the result to mergings.
    """
    #mergings = linkage(samples, method='complete')
    mergings = linkage(atoms, method='complete')

    dendrogram(mergings,
               leaf_rotation=90,
               leaf_font_size=6,
               )
    
    plt.show()
    
def tSNE(atoms):
    # Importing Modules
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt
    
    # Loading dataset
    #iris_df = datasets.load_iris()
    
    # Defining Model
    model = TSNE(learning_rate=100)
    
    # Fitting Model
    transformed = model.fit_transform(atoms)
    print transformed
    
    # Plotting 2d t-Sne
    x_axis = transformed[:, 0]
    y_axis = transformed[:, 1]
    
    plt.scatter(x_axis, y_axis)#, c=iris_df.target)
    plt.show()
#---------------------------------------------------------------------------
# Purpose: Plotting scatter plot
def visualize_coordinates(x,y,z):
    fig = plt.figure() # define figure object
    ax = Axes3D(fig)
    ax.plot(x, y, z) # scatter plot
    plt.show() # show interactive image
#---------------------------------------------------------------------------
# Main
#---------------------------------------------------------------------------
if __name__ == "__main__":
    
    ### Getting files information and loading data into lists
    psfFile = "/home/jordancaraballo/Documents/Research/Cellulose/Tetradecane/Tetradecane100/Tetradecane100.psf" # name of the psf file
    psfList = parsePSF_toList(psfFile) # takes the list of psf lines
    
    pdbFile = "/home/jordancaraballo/Documents/Research/Cellulose/Tetradecane/Tetradecane100/Tetradecane100_lastFrame.pdb"
    pdbList = parsePDB_toList(pdbFile)

    ### Starting to define spacer arms information and coordinates
    pabaAtoms = getAtomsMol(psfList, 'PAB') # takes paba atoms from psflist
    baseAtoms  = getAtoms(psfList, 'C15', tcolumn=5) # takes C15 (base carbon) atom from each molecule
    upperAtoms = getAtoms(psfList, 'C60', tcolumn=5) # takes C06 (upper carbon) atom from each molecule
    
    ### Starting to define graph and shortests path
    numBonds  = (subprocess.check_output("cat " + psfFile + "| grep NBOND | cut -d \" \" -f4", shell=True))[:-1]
    bondConnects = getBonds(psfList, numBonds) # get bond connections list
    
    G = nx.Graph()
    G.add_edges_from(bondConnects)
    Gsub = G.subgraph(pabaAtoms)
    #nx.draw(Gsub)

    ### Starting to get atom indices and coordinates from the spacer arm    
    comp = list(Gsub.subgraph(c) for c in nx.connected_components(Gsub))
    graphIndexes = getSpacerArmGraph(comp,baseAtoms,upperAtoms) #list of atom indexes
    print (len(graphIndexes[0]))
    atomCoordinates = getPABACoordinates(graphIndexes, psfList, pdbList) # numpy Array

    reshapedCoordinates = list()
    for molecule in atomCoordinates:
        reshape = molecule.reshape(3,12)
        reshapedCoordinates.append(np.concatenate([reshape[0], reshape[1], reshape[2]], axis=0))
    reshapedCoordinates = np.array(reshapedCoordinates)
    
    kmeans(reshapedCoordinates)
    #hierarchical(reshapedCoordinates)
    tSNE(reshapedCoordinates)
    ### For each molecule in the system
    #coordinates_rotated = list() # will be a list of numpy arrays
    #for molecule in atomCoordinates:
        #plt.plot(molecule[:,0], molecule[:,1], 'ro')
    #    visualize_coordinates(molecule[:,0],molecule[:,1], molecule[:,2])
    #    plt.show()
        
        #translated     = translate_coordinates_lowestAtom(molecule)  # translate
        #rotated_RyRx   = rotate_matrix_RyRxHighest(translated)       # rotate RyRx
        #rotated_RyRxRz = rotate_matrix_RzMiddle(rotated_RyRx)        # rotate Rz
        #coordinates_rotated.append(rotated_RyRxRz)                   # append modified coordinates
        #getVarSpiral_A(rotated_RyRxRz)                               # get A from spiral equation
        #getVarSpiral_R(rotated_RyRxRz)                               # get R from spiral equation
    
        ### Spiral equation variables
        # Parametric Equation: E = (R cos P t + x0 + (alpha)(t), R sen P t + y0 + (beta)(t), (A)(t) + z0 + (lamda)(t))
        
    """
    # Pasos necesarios:
    1. Definir funcion F([R,x0,alpha,y0,betha,z0,lamda,P])
    2. min t = (Ex(t)-xi)^2 + (Ey(t)-yi)^2 + (Ez(t)-zi)^2 ===  con scipy
    3. Resolver (d (min t)/dt) = 0 con algun metodo numero (Newton maybe)
    """
    ### Solving for the Spiral equation
    
    
    
    
    
    
    
    
    
