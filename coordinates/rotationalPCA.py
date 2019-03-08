"""
Describir proposito luego
"""
import numpy as np # include numpy
import math
from numpy import genfromtxt # read from csv file

# export file coordinates into numpy array
coors = genfromtxt('/home/jordancaraballo/Documents/Research/Cellulose/functionalizedCellulose/data/TetradecanePCA/RotationalPCA/Tetradecane01_LastFramePABACoordinates.txt', delimiter=',')

# split arrays into n arrays (Tetradecane005 = 3 spacer arms = 3 arrays)
coors_splitted = np.split(coors, 3)

# indices de los atomos que forman el triangulo = [C6: 26, C11: 38, C15: 20]
# llevar coordenadas al origen con la resta de C15 (atomo inferior)
coors_updated_origin = [coors_splitted[0] - coors_splitted[0][20], coors_splitted[1] - coors_splitted[1][20], coors_splitted[2] - coors_splitted[2][20]]

# formar triangulo entre los puntos escogidos
V1_C6,V2_C6,V3_C6    = coors_updated_origin[0][26] # matriz RyRxC6
print ("Updated C6 (V1,V2,V3): ", V1_C6,V2_C6,V3_C6)

# Solving for rotational matrix
# a = cos x, b = sin x, c = cos x, d = sin x
a = math.sqrt(V3_C6**2 / (V2_C6**2 + V3_C6**2))
b = math.sqrt(1 - V3_C6**2 / (V2_C6**2 + V3_C6**2))
c = math.sqrt(1 - V1_C6**2 / ((V3_C6 * a + V2_C6 * b)**2 + V1_C6**2))
d = math.sqrt(V1_C6**2 / ((V3_C6 * a + V2_C6 * b)**2 + V1_C6**2))
print ("Updated C6 (a: ", a , "b: ", b, "c: ", c, "d: ", d, ")")

# definiendo matrices de rotacion generales
RyRxC6_matrix = np.array([[c, d * b, d * a], [0, a, -b], [-d, b * c, a * c]])
coors_updated_RyRx = coors_updated_origin[0].dot(RyRxC6_matrix)

# definiendo la ultima rotacion RzC11
V1_C11,V2_C11,V3_C11 = coors_updated_RyRx[38] # Por el momento solo se esta tomando una molecula

# e = cos z, f = sin z
e = math.sqrt(V2_C11**2 / (V1_C11**2 + V2_C11**2))
f = math.sqrt(1 - V2_C11**2 / (V1_C11**2 + V2_C11**2))
RzC11_matrix = np.array([[e, -f, 0], [f, e, 0], [0, 0, 1]])
coors_updated_Rz = coors_updated_RyRx.dot(RzC11_matrix)

print (coors_updated_Rz)


# Plotting
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')
#x.scatter3D(coors_updated_origin[0][:,0], coors_updated_origin[0][:,1], coors_updated_origin[0][:,2])
#ax.scatter3D(coors_splitted[1][:,0], coors_splitted[1][:,1], coors_splitted[1][:,2])
#ax.scatter3D(coors_updated_RyRx[:,0], coors_updated_RyRx[:,1], coors_updated_RyRx[:,2])
ax.scatter3D(coors_updated_Rz[:,0], coors_updated_Rz[:,1], coors_updated_Rz[:,2])

plt.show()































"""
# Rotational Matrix
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

ax = plt.axes(projection='3d')

# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)
#ax.plot3D(xline, yline, zline, 'gray')
#ax.plot3D(coors[:,[0]], coors[:,[1]], coors[:,[2]], 'gray')

# Data for three-dimensional scattered points
#zdata = 15 * np.random.random(100)
#xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
#ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
#ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');

ax.scatter3D(coors[:,0], coors[:,1], coors[:,2])#, c=coors[:,2], cmap='Greens');
ax.plot3D(coors[:,0], coors[:,1], coors[:,2])#, c=coors[:,2], cmap='Greens');

plt.show()
"""
