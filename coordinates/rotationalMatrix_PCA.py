#############################################################################
# Purpose: Produce 3-dimensional rotational matrix to rotate molecule
# coordinates into the same plane (x,y,z). This will ease PCA calculations.
# Output are rotational matrixs, updated coordinates, and a PCA.
# Author: Jordan A Caraballo Vega
# Collaborator: Jose O Sotero Esteva, Pablo Negron
# Department of Mathematics, University of Puerto Rico at Humacao
#############################################################################
# Imports
import math                      # for math functions
import numpy as np               # for numpy arrays and calculations
from numpy import genfromtxt     # read from csv file
import matplotlib.pyplot as plt  # for visualizations
from mpl_toolkits import mplot3d # for visualizations
from sklearn.decomposition import PCA # for PCA
#---------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------
# Purpose: Plotting scatter plot
def visualize_coordinates(x,y,z):
    fig = plt.figure() # define figure object
    ax  = plt.axes(projection='3d') # define three dimensional property
    ax.scatter3D(x, y, z) # scatter plot
    plt.show() # show interactive image
#---------------------------------------------------------------------------
# Purpose: translate coordinates to the origin based on lowest pABA atom
def translate_coordinates_lowestAtom(molecule_array, index_lowest):
    return molecule_array - molecule_array[index_lowest]
#---------------------------------------------------------------------------
# Purpose: rotate coordinates over RyRx axis
def rotate_matrix_RyRxHighest(molecule, index_highest):
    V1,V2,V3 = molecule[index_highest] # coordinates defined as vector representations

    # Solving for rotational matrix: a = cos x, b = sin x, c = cos x, d = sin x
    a = math.sqrt(V3**2 / (V2**2 + V3**2))
    b = math.sqrt(1 - V3**2 / (V2**2 + V3**2))
    c = math.sqrt(1 - V1**2 / ((V3 * a + V2 * b)**2 + V1**2))
    d = math.sqrt(V1**2 / ((V3 * a + V2 * b)**2 + V1**2))

    return molecule.dot(np.array([[c, d * b, d * a], [0, a, -b], [-d, b * c, a * c]])) # return rotated molecule matrix
#---------------------------------------------------------------------------
# Purpose: rotate coordinates over Rz axis
def rotate_matrix_RzMiddle(molecule, index_middle):
    V1,V2,V3 = molecule[index_middle] # coordinates defined as vector representations

    # Solving for rotational matrix: e = cos z, f = sin z
    e = math.sqrt(V2**2 / (V1**2 + V2**2))
    f = math.sqrt(1 - V2**2 / (V1**2 + V2**2))

    return molecule.dot(np.array([[e, -f, 0], [f, e, 0], [0, 0, 1]])) # return rotated molecule matrix

def plotting_pca(molecules, dim):
    # concatenating arrays
    system_array = molecules[0]
    for molecule in molecules[1:]:
        system_array = np.concatenate((system_array,molecule),axis=0)

    pca = PCA(n_components=dim)
    projected = pca.fit_transform(system_array)
    #projected = pca.inverse_transform(projected)

    if dim == 2:
        plt.scatter(projected[:, 0], projected[:, 1], alpha=0.8)
    if dim == 3:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3D(projected[:, 0], projected[:, 1], projected[:, 2], alpha=0.5, cmap=plt.cm.get_cmap('Accent', 10), edgecolor='none',)
    plt.show()
#---------------------------------------------------------------------------
# Main
#---------------------------------------------------------------------------
if __name__ == "__main__":

    # export file coordinates into numpy array
    file_input = genfromtxt('/home/jordancaraballo/Documents/Research/Cellulose/functionalizedCellulose/'+ \
                       'data/TetradecanePCA/RotationalPCA/Tetradecane005_LastFramePABACoordinates.txt', delimiter=',')
    n_atoms = 52 # number of atoms per spacer arm
    # list index of atoms forming triangle for rotation = [C6: 26, C11: 38, C15: 20]
    # 005, 01 = [C6: 26, C11: 38, C15: 20]
    # 02 = [C6: 29, C11: 38, C15: 33], 43
    triangle_atoms_indexes_tetra = {'lowest': 26, 'middle': 38, 'highest': 20} # dictionary with indexes of triangle atoms

    # split arrays into n arrays (Tetradecane005 = 3 spacer arms = 3 arrays), list of numpy arrays
    coordinates = np.split(file_input, len(file_input) / n_atoms) # (# atoms) / (# atoms spacer arm) = # of spacer arms
    print ("Molecules: ", len(coordinates))
    coordinates_rotated = list() # store rotated molecules here, list of numpy arrays

    # Iterate through list of numpy arrays
    for molecule in coordinates:
        translated     = translate_coordinates_lowestAtom(molecule, triangle_atoms_indexes_tetra['lowest'])  # translate
        rotated_RyRx   = rotate_matrix_RyRxHighest(translated,      triangle_atoms_indexes_tetra['highest']) # rotate RyRx
        rotated_RyRxRz = rotate_matrix_RzMiddle(rotated_RyRx,       triangle_atoms_indexes_tetra['middle'])  # rotate Rz
        coordinates_rotated.append(rotated_RyRxRz)

    # Testing
    #visualize_coordinates(coordinates[1][:,0], coordinates[1][:,1], coordinates[1][:,2])
    #visualize_coordinates(coordinates_rotated[1][:,0], coordinates_rotated[1][:,1], coordinates_rotated[1][:,2])
    plotting_pca(coordinates,3)
    plotting_pca(coordinates_rotated,3)
