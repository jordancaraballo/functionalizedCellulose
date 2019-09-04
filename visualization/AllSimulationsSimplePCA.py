
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
