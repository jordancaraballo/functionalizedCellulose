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
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math                                        # import for distance calculation
import numpy as np                                 # include numpy

from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
#----------------------------------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------------------------------
def parsePSF_toList(psfFile): # input: string with psf file name, output: list of file lines splitted
    print "Parsing PSF: " + psfFile
    return [line.strip().split() for line in open(psfFile,'r').readlines() ][6:] # list with psf file lines
#----------------------------------------------------------------------------------------------------------------
def parsePDB_toList(pdbFile): # input: string with pdb file name, output: list of file lines splitted
    print "Parsing PDB: " + pdbFile
    return [line.strip().split() for line in open(pdbFile,'r').readlines() ][1:] # list with pdb file lines
#----------------------------------------------------------------------------------------------------------------
def parseDCD_toList(dcdFile): # input: string with dcd file name, output: list of dcds copies
    print "Parsing DCD..."
    return [DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile), DCDReader(dcdFile)] # list with dcd files
#----------------------------------------------------------------------------------------------------------------
# combination of pdb and psf file to parse paba coordinates
# input: dcd, pdb list, psf list, residue to match. output: list of atoms
def getAtoms(psf_list, residue, tcolumn=5, mcolumn=3, mol="PAB"):
    #print "Extracting Atoms for residue " + residue + " In mol " + mol + "..."
    atomsInfo = list() # list to store atoms
    for atom in psf_list: # for range in the number of atoms
        try:
            if atom[tcolumn][:len(residue)] == residue and atom[mcolumn][:len(mol)] == mol: # if atom is PABA go for it
                atomsInfo.append(int(atom[0])) # append atom line to list
        except IndexError:
            pass
    print "Total Atoms: " + str(len(atomsInfo)) + " for residue " + residue + " In mol " + mol
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
    print "Total Atoms: " + str(len(atomsInfo)) + " from mol " + mol
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

def getNitrogenCoordinates(indexesList, psfList, pdbList):
    atomsCoordinates = list()
    for nit in indexesList:
        atomsCoordinates.append(float(pdbList[nit-1][8]))
    return atomsCoordinates
    
def reshapeCoordinates(atomCoordinates,x=3,y=12,dim='3D'):
    reshapedCoordinates = list()
    for molecule in atomCoordinates:
        reshape = molecule.reshape(x,y)
        if dim == '2D':
            reshapedCoordinates.append(np.concatenate([reshape[0], reshape[1]], axis=0))
        else:
            reshapedCoordinates.append(np.concatenate([reshape[0], reshape[1], reshape[2]], axis=0))
    return np.array(reshapedCoordinates)
    
def getDistancesAngles(atomCoordinates):
    newDistances = list()
    for molecule in atomCoordinates:
        d = np.diff(molecule, axis=0)
        segdists = np.sqrt((d ** 2).sum(axis=1)) # numkpy array with distances between adjacent points
        print "Segdists Size: ", segdists.shape
        
        angles = list()
        for atom in range(len(molecule)-2):
            f = molecule[atom+1] - molecule[atom]   # b-a 
            e = molecule[atom+1] - molecule[atom+2] # b-c 
            abNorm = f / norm(f)
            bcNorm = e / norm(e)
            angle = np.arccos(abNorm[0] * bcNorm[0] + abNorm[1] * bcNorm[1] + abNorm[2] * bcNorm[2])*180.0 / np.pi
            angles.append(angle)
        angles.append(angles[-1])
        newDistances.append(np.column_stack((segdists, np.array(angles))))
    return np.array(newDistances)
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
def rotate_matrix_full(molecules, low=11, mid=5, high=0, index_middle=5):
    coordinates_rotated = list()
    for mol in molecules:
        translated     = translate_coordinates_lowestAtom(mol, low)  # translate
        rotated_RyRx   = rotate_matrix_RyRxHighest(translated, high) # rotate RyRx
        rotated_RyRxRz = rotate_matrix_RzMiddle(rotated_RyRx,  mid)  # rotate Rz
        coordinates_rotated.append(rotated_RyRxRz)
    return np.array(coordinates_rotated)
#---------------------------------------------------------------------------
# Purpose: Find A of the spiral equation = A: height of the upper atom
def getVarSpiral_A(molecule, index_highest=-1, index_base=0):
    return molecule[index_highest][2] - molecule[index_base][2]
#---------------------------------------------------------------------------
def kmeans(atoms, colors):
    model = KMeans(n_clusters=3)
    model.fit(atoms)

    #all_predictions = model.predict(atoms)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #scatter = ax.scatter(atoms[:,0],atoms[:,1], c=all_predictions,s=50)
    #scatter = ax.scatter(atoms[:,0],atoms[:,1], c=colors)#,s=50)
    ax.scatter(atoms[:,0],atoms[:,1], c=colors)#,s=50)
    ax.set_title('K-Means Clustering')
    #plt.colorbar(scatter)
    #print(all_predictions)
    plt.show()
    ax.savefig('kmeans.png')

def kmeansNew(atoms):
    model = KMeans(n_clusters=3)
    model.fit(atoms)
    
    # Prediction on the entire data
    all_predictions = model.predict(atoms)
    print atoms[:72,0].shape
    print all_predictions.shape
    
    #Plot the clusters obtained using k means
    fig = plt.figure()
    ax = fig.add_subplot(111)
    scatter = ax.scatter(atoms[:,0],atoms[:,1], c=all_predictions,s=50)
    ax.set_title('K-Means Clustering')
    plt.colorbar(scatter)
    print(all_predictions)
    plt.savefig('kmeansnew.png')

    
def hierarchical(atoms):
    #mergings = linkage(samples, method='complete')
    mergings = linkage(atoms, method='complete')
    dendrogram(mergings,
               leaf_rotation=90,
               leaf_font_size=6,
               )    
    plt.show()
  
def tSNE(atoms, colors, learningRate=100):
    model = TSNE(learning_rate=learningRate)
    transformed = model.fit_transform(atoms)
    x_axis = transformed[:, 0]
    y_axis = transformed[:, 1]
    #plt.scatter(x_axis, y_axis,c=atoms[:,2])#, c=iris_df.target)
    plt.scatter(x_axis, y_axis,c=colors)#, c=iris_df.target)
    plt.title('tSNE Clustering')
    plt.show()
    plt.savefig('tsne.png')

def dbscan(atoms, colors):
    dbscan = DBSCAN()
    dbscan.fit(atoms)
    pca = PCA(n_components=2).fit(atoms)
    pca_2d = pca.transform(atoms)
    # Plotting 2d t-Sne
    x_axis = pca_2d[:, 0]
    y_axis = pca_2d[:, 1]
    plt.scatter(x_axis, y_axis, c=colors)#, c=iris_df.target)
    plt.title('dbscan Clustering')
    plt.show()
    #plt.savefig('dbscan.png')

def pca2d(atoms, colors):
    from sklearn.decomposition import PCA
    from matplotlib.lines import Line2D
    
    ### Calculate eigenvectors
    pca = PCA(n_components=2)
    X_reduced =pca.fit_transform(atoms)
    plt.title('PCA2D Clustering')
    plt.scatter(X_reduced[:, 0],X_reduced[:, 1], c=colors, alpha=0.8, edgecolor='none')
  
    base_colors = ['blue','green','red','cyan','magenta','yellow','black','indigo','burlywood','darksalmon','darkviolet']
    legend_elements = [ Line2D([0], [0], marker='o', color=base_colors[0],  label='5%',   markerfacecolor=base_colors[0], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[1],  label='10%',  markerfacecolor=base_colors[1], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[2],  label='20%',  markerfacecolor=base_colors[2], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[3],  label='30%',  markerfacecolor=base_colors[3], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[4],  label='40%',  markerfacecolor=base_colors[4], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[5],  label='50%',  markerfacecolor=base_colors[5], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[6],  label='60%',  markerfacecolor=base_colors[6], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[7],  label='70%',  markerfacecolor=base_colors[7], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[8],  label='80%',  markerfacecolor=base_colors[8], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[9],  label='90%',  markerfacecolor=base_colors[9], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[10], label='100%', markerfacecolor=base_colors[10], markersize=10),
    ]
    plt.legend(loc='upper left', handles=legend_elements, frameon=False) 
    plt.show()
    plt.savefig('pca2d.png')
    
def pca2d_OnClick(atoms, colors):
    from sklearn.decomposition import PCA
    from matplotlib.lines import Line2D

    ### Calculate eigenvectors
    pca = PCA(n_components=2)
    X_reduced = pca.fit_transform(atoms)
        
    # draw a scatter plot of the generated values
    fig = plt.figure(figsize=(20, 16))
    ax = fig.add_subplot(111)

    # legend elements
    base_colors = ['blue','green','red','cyan','magenta','yellow','black','indigo','burlywood','darksalmon','darkviolet']
    legend_elements = [ Line2D([0], [0], marker='o', color=base_colors[0],  label='5%',   markerfacecolor=base_colors[0], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[1],  label='10%',  markerfacecolor=base_colors[1], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[2],  label='20%',  markerfacecolor=base_colors[2], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[3],  label='30%',  markerfacecolor=base_colors[3], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[4],  label='40%',  markerfacecolor=base_colors[4], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[5],  label='50%',  markerfacecolor=base_colors[5], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[6],  label='60%',  markerfacecolor=base_colors[6], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[7],  label='70%',  markerfacecolor=base_colors[7], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[8],  label='80%',  markerfacecolor=base_colors[8], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[9],  label='90%',  markerfacecolor=base_colors[9], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[10], label='100%', markerfacecolor=base_colors[10], markersize=10),
    ]
    plt.legend(loc='upper left', handles=legend_elements, frameon=False) 
    plt.title('Para-aminobenzamidine Morphology - Principal Components Analysis 14-atoms')
    
    # extract the scatterplot drawing in a separate function so we ca re-use the code
    ax.scatter( X_reduced[:, 0], X_reduced[:, 1], c=colors, alpha=0.8, edgecolor='none', s=50, picker=True)

    # define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
    def onpick(event):
        # Create embedded figure            
        figi = plt.figure()
        for subplotnum, dataind in enumerate(event.ind):
            ax = Axes3D(figi)
            data = atoms[dataind].reshape(12,3)
            ax.plot(data[:,0],data[:,1],data[:,2],c=colors[dataind])
            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Y Axis')
            ax.set_title('Spacer Arm Coordinates')
        figi.show()

    # connect the click handler function to the scatterplot
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show()
    fig.savefig('pca2d-interactive.png')

def pca2d_OnClickFit(atoms, colors):
    from sklearn.decomposition import PCA
    from matplotlib.lines import Line2D
    #from numpy import arange, cos, linspace, pi, sin, random
    from scipy.interpolate import splprep, splev

    ### Calculate eigenvectors
    pca = PCA(n_components=2)
    X_reduced = pca.fit_transform(atoms)
        
    # draw a scatter plot of the generated values
    fig = plt.figure(figsize=(20, 16))
    ax = fig.add_subplot(111)

    # legend elements
    base_colors = ['blue','green','red','cyan','magenta','yellow','black','indigo','burlywood','darksalmon','darkviolet']
    legend_elements = [ Line2D([0], [0], marker='o', color=base_colors[0],  label='5%',   markerfacecolor=base_colors[0], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[1],  label='10%',  markerfacecolor=base_colors[1], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[2],  label='20%',  markerfacecolor=base_colors[2], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[3],  label='30%',  markerfacecolor=base_colors[3], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[4],  label='40%',  markerfacecolor=base_colors[4], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[5],  label='50%',  markerfacecolor=base_colors[5], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[6],  label='60%',  markerfacecolor=base_colors[6], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[7],  label='70%',  markerfacecolor=base_colors[7], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[8],  label='80%',  markerfacecolor=base_colors[8], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[9],  label='90%',  markerfacecolor=base_colors[9], markersize=10),
                        Line2D([0], [0], marker='o', color=base_colors[10], label='100%', markerfacecolor=base_colors[10], markersize=10),
    ]
    plt.legend(loc='upper left', handles=legend_elements, frameon=False) 
    plt.title('Para-aminobenzamidine Morphology - Principal Components Analysis 5-atoms')
    
    # extract the scatterplot drawing in a separate function so we ca re-use the code
    ax.scatter( X_reduced[:, 0], X_reduced[:, 1], c=colors, alpha=0.8, edgecolor='none', s=50, picker=True)

    # define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
    def onpick(event):
        # Create embedded figure 
        # spline parameters
    
        s=3.0 # smoothness parameter
        k=2 # spline order
        
        figi = plt.figure()
        for subplotnum, dataind in enumerate(event.ind):
            ax = Axes3D(figi)
            data = atoms[dataind].reshape(12,3)
            tckp,u = splprep([data[:,0],data[:,1],data[:,2]],s=s,k=k,nest=-1)   # find the knot points 
            xnew,ynew,znew = splev(np.linspace(0,1,400),tckp)
            print xnew,ynew,znew
            ax.plot(data[:,0],data[:,1],data[:,2],'o', c=colors[dataind])
            ax.plot(xnew,ynew,znew,'r-',label='fit',c='k')
            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Y Axis')
            ax.set_title('Spacer Arm Coordinates')
        figi.show()

    # connect the click handler function to the scatterplot
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show()
    fig.savefig('pca2d-interactive.png')
    
def pca3d(atoms, colors):
    from sklearn.decomposition import PCA
    fig = plt.figure(1, figsize=(8, 6))
    ax = Axes3D(fig, elev=-150, azim=110)
    pca = PCA(n_components=3)
    X_reduced =pca.fit_transform(atoms)
    ax.scatter3D(X_reduced[:, 0],X_reduced[:, 1],X_reduced[:, 2], c=colors)
    plt.title('PCA3D Clustering')
    plt.show()
    #plt.savefig('pca3d.png')
#---------------------------------------------------------------------------
# Purpose: Plotting scatter plot
def visualize_coordinates(x,y,z):
    fig = plt.figure() # define figure object
    ax = Axes3D(fig)
    ax.plot(x, y, z) # scatter plot
    plt.show() # show interactive image
#---------------------------------------------------------------------------
def frequencyChart(data, dMin, dMax, colors, system):
    # Define figures for the plot
    yValues = [5.00, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.00]

    fig = plt.figure()
    fig.suptitle('Distribution of PABA Heights',fontsize=18)
    ax  = fig.add_subplot(111)
    fig.subplots_adjust(top=0.90)
    ax.set_title('14-atoms spacer arm')
    ax.set_xlabel('PABA Height (A)')
    ax.set_ylabel('PABA Occupancy')
    # Filter and calculate frequencies
    for i in range(data.shape[0]):
    
        # Calculate frequencies
        h = np.histogram(data[i],20,(dMin,dMax),density=True)
        n = len(h[1][:-1])
        x = h[1][:-1]
        y = np.ones(n)*yValues[i]
        sizes=(h[0]*5000)
        ax.scatter(x,y,s=sizes,c=colors[i])
        ax.set_xlim([15, 30])
    
    plt.show()
    fig.savefig('frequency'+system+'.png')
