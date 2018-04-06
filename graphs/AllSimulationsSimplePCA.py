'''
import pandas as pd
import plotly.plotly as py
py.sign_in('jordancaraballo', 's0DmyorY934VFkO6jBdA')
from plotly.graph_objs import *
import plotly.tools as tls
from matplotlib.mlab import PCA
import numpy
import matplotlib.pyplot as plt

data = pd.read_csv( filepath_or_buffer='Tetradecane05PCA.txt', header=None, sep=',')

#print df

#df = pd.read_csv(
#    filepath_or_buffer='https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data',
#    header=None,
#    sep=',')

def doPCA(data):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=52)
    pca.fit(data)
    return pca

pca = doPCA(data)
print pca.explained_variance_ratio_
first_pc = pca.components_[0]
second_pc = pca.components_[1]

transformed_data = pca.transform(data)
for ii, jj in zip(transformed_data, data):
    plt.scatter( first_pc[0]*ii[0], first_pc[1]**ii[0], color="r")
    plt.scatter( second_pc[0]*ii[1], second_pc[1]**ii[1], color="c")
    #plt.scatter( jj[0], jj[1], color="b")

plt.xlabel("yuu")
plt.ylabel("yoo")
plt.show()

import pandas as pd
import plotly.plotly as py
py.sign_in('jordancaraballo', 's0DmyorY934VFkO6jBdA')
from plotly.graph_objs import *
import plotly.tools as tls
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(4294967295) # random seed for consistency

# A reader pointed out that Python 2.7 would raise a
# "ValueError: object of too small depth for desired array".
# This can be avoided by choosing a smaller random seed, e.g. 1
# or by completely omitting this line, since I just used the random seed for
# consistency.

mu_vec1 = np.array([0,0,0])
cov_mat1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
class1_sample = np.random.multivariate_normal(mu_vec1, cov_mat1, 20).T
assert class1_sample.shape == (3,20), "The matrix has not the dimensions 3x20"

mu_vec2 = np.array([1,1,1])
cov_mat2 = np.array([[1,0,0],[0,1,0],[0,0,1]])
class2_sample = np.random.multivariate_normal(mu_vec2, cov_mat2, 20).T
assert class2_sample.shape == (3,20), "The matrix has not the dimensions 3x20"

#data = pd.read_csv( filepath_or_buffer='Tetradecane005PCA.txt', header=None, sep=',')
data = np.loadtxt('Tetradecane005PCA.txt',delimiter=',')
#print data
#print class1_sample

class1_sample = data
#print class1_sample


#print " NOOOJOJO "
#print class1_sample[0,:]
#print class1_sample[1,:]



from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
plt.rcParams['legend.fontsize'] = 10
#ax.plot(class1_sample[0,:], class1_sample[1,:], class1_sample[2,:], 'o', markersize=8, color='blue', alpha=0.5, label='class1')
#ax.plot(class2_sample[0,:], class2_sample[1,:], class2_sample[2,:], '^', markersize=8, alpha=0.5, color='red', label='class2')
#ax.plot(class1_sample[0,:], class1_sample[1,:], class1_sample[2,:], 'o', markersize=8, color='blue', alpha=0.5, label='class1')

ax.plot(class1_sample[0,:], class1_sample[0,:],class1_sample[0,:],'o', markersize=8, color='blue', alpha=0.5, label='class1')

plt.title('Samples for class 1 and class 2')
ax.legend(loc='upper right')

plt.show()
'''

from sklearn.decomposition import PCA as sklearnPCA
import matplotlib.pyplot as plt
import numpy as np

data_1 = np.loadtxt('Tetradecane005PCA2.txt',delimiter=',')
data_2 = np.loadtxt('Tetradecane01PCA2.txt',delimiter=',')
data_3 = np.loadtxt('Tetradecane02PCA2.txt',delimiter=',')
data_4 = np.loadtxt('Tetradecane03PCA2.txt',delimiter=',')
data_5 = np.loadtxt('Tetradecane04PCA2.txt',delimiter=',')
data_6 = np.loadtxt('Tetradecane05PCA2.txt',delimiter=',')
data_7 = np.loadtxt('Tetradecane06PCA2.txt',delimiter=',')
data_8 = np.loadtxt('Tetradecane07PCA2.txt',delimiter=',')
data_9 = np.loadtxt('Tetradecane08PCA2.txt',delimiter=',')
#data_10 = np.loadtxt('Tetradecane09PCA2.txt',delimiter=',')
data_11 = np.loadtxt('Tetradecane100PCA2.txt',delimiter=',')

#all_samples = np.concatenate((data_1, data_2), axis=1)
#print all_samples.shape

sklearn_pca = sklearnPCA(n_components=2)

sklearn_transf1 = sklearn_pca.fit_transform(data_1.T)
sklearn_transf2 = sklearn_pca.fit_transform(data_2.T)
sklearn_transf3 = sklearn_pca.fit_transform(data_3.T)
sklearn_transf4 = sklearn_pca.fit_transform(data_4.T)
sklearn_transf5 = sklearn_pca.fit_transform(data_5.T)
sklearn_transf6 = sklearn_pca.fit_transform(data_6.T)
sklearn_transf7 = sklearn_pca.fit_transform(data_7.T)
sklearn_transf8 = sklearn_pca.fit_transform(data_8.T)
sklearn_transf9 = sklearn_pca.fit_transform(data_9.T)
#sklearn_transf10 = sklearn_pca.fit_transform(data_10.T)
sklearn_transf11 = sklearn_pca.fit_transform(data_11.T)

plt.plot(sklearn_transf1[0:156,0], sklearn_transf1[0:156,1], 'o', markersize=7, color='blue', alpha=0.5, label='5%')
plt.plot(sklearn_transf2[0:156,0], sklearn_transf2[0:156,1], 'o', markersize=7, color='red', alpha=0.5, label='10%')
plt.plot(sklearn_transf3[0:156,0], sklearn_transf3[0:156,1], 'o', markersize=7, color='green', alpha=0.5, label='20%')
plt.plot(sklearn_transf4[0:156,0], sklearn_transf4[0:156,1], 'o', markersize=7, color='orange', alpha=0.5, label='30%')
plt.plot(sklearn_transf5[0:156,0], sklearn_transf5[0:156,1], 'o', markersize=7, color='black', alpha=0.5, label='40%')
plt.plot(sklearn_transf6[0:156,0], sklearn_transf6[0:156,1], 'o', markersize=7, color='purple', alpha=0.5, label='50%')
plt.plot(sklearn_transf7[0:156,0], sklearn_transf7[0:156,1], 'o', markersize=7, color='grey', alpha=0.5, label='60%')
plt.plot(sklearn_transf8[0:156,0], sklearn_transf8[0:156,1], 'o', markersize=7, color='yellow', alpha=0.5, label='70%')
plt.plot(sklearn_transf9[0:156,0], sklearn_transf9[0:156,1], 'o', markersize=7, color='pink', alpha=0.5, label='80%')
#plt.plot(sklearn_transf10[0:156,0], sklearn_transf10[0:156,1], 'o', markersize=7, color='cyan', alpha=0.5, label='90%')
plt.plot(sklearn_transf11[0:156,0], sklearn_transf11[0:156,1], 'o', markersize=7, color='brown', alpha=0.5, label='100%')


plt.xlabel('x_values')
plt.ylabel('y_values')
#plt.xlim([-4,4])
#plt.ylim([-4,4])
plt.legend()
plt.title('Transformed samples with class labels from matplotlib.mlab.PCA()')

plt.show()
