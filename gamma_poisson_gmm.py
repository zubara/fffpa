# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:29:26 2016

@author: zubarei1
"""

#import numpy as np
#from time import time
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.ticker import NullFormatter
#from scipy.io import loadmat
#from bayespy.nodes import Poisson
#from xlrd import open_workbook
#import os
#os.chdir('/m/home/home6/62/zubarei1/data/Desktop/elinas_APs')
#a = loadmat('abfdata_08-Jan-2017.mat')
#X = a['X']#[:,22:22+18]
#spikes = a['spikes']
#sweeps = a['X_2sw']
##cells = a['cellstats']

X-=X.mean(0)
X/=X.std(0)
from sklearn.manifold import TSNE, MDS, SpectralEmbedding, LocallyLinearEmbedding
#mlab = loadmat('raw_man_lab.mat')
#labelsm = mlab['man_labels']

from bayespy.nodes import Dirichlet, Categorical
from bayespy.nodes import Gaussian, Wishart, GaussianARD, Gamma
from bayespy.nodes import Mixture
    #from bayespy.nodes import GaussianARD, Dirichlet, GaussianWishart, Categorical
from bayespy.inference import VB
"""GMM MODEL"""
plt.close('all')
sc_alpha = 1e-3
sc_m0 = 1e-3
sc_gam_rate = 1e-4
sc_gam_shape = 1e-5

n, d = X.shape
k = 10
m0 = GaussianARD(0,sc_m0,shape=(d,), plates=(k,), name='m0')
lam = Gamma(sc_gam_shape,sc_gam_rate,plates=(k,d,), name='lam')

#g = GaussianARD(m0,lam)

a0 = sc_alpha*np.ones(k)/k
alpha = Dirichlet(a0, name='alpha')
Z = Categorical(alpha, plates=(n,1), name='Z')
Y = Mixture(Z, GaussianARD, m0, lam, name='Y',cluster_plate=-2)

Z.initialize_from_random()
Q = VB(Y, m0, lam, Z, alpha)
Y.observe(X)
Q.update(repeat=1000)
means = Q['m0'].get_moments()[0]
precs = Q['lam'].get_moments()[0]
alphas = Q['alpha'].get_moments()[0]
Zs = Q['Z'].get_moments()[0]
print(np.exp(alphas))
act_comp = np.where(np.exp(alphas)>1e-6)[0]
Zc = np.squeeze(Zs)
Zc = Zc[...,act_comp]
zmax = np.argmax(Zc,-1)

f,ax = plt.subplots(2,3)
n_neighbors = 5
n_components = 2


colors = np.array(['r', 'g', 'b', 'y', 'c', 'm', 'k'])#, [0.5,0.5,0.5], [0.75, 0.75,0]])
t0 = time()
mds = MDS(n_components, max_iter=1000, n_init=1)
Y = mds.fit_transform(X)
t1 = time()
print("MDS: %.2g sec" % (t1 - t0))
#ax = fig.add_subplot(258)
#a[].figure()
ax[0,0].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
ax[0,0].set_title("MDS (%.2g sec)" % (t1 - t0))

t0 = time()
se = SpectralEmbedding(n_components=n_components, n_neighbors=n_neighbors)
Y = se.fit_transform(X)
t1 = time()
print("SpectralEmbedding: %.2g sec" % (t1 - t0))
ax[0,1].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
ax[0,1].set_title("SpectralEmbedding (%.2g sec)" % (t1 - t0))

Xp = X - X.min(0)
t0 = time()
tsne = TSNE(n_components=n_components, init='pca', random_state=42)
Y = tsne.fit_transform(Xp)
t1 = time()
print("t-SNE: %.2g sec" % (t1 - t0))
ax[0,2].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
ax[0,2].set_title("t-SNE (%.2g sec)" % (t1 - t0))

methods = ['standard', 'modified', 'ltsa']
labels = ['LLE', 'Modified LLE', 'LTSA']

for i, method in enumerate(methods):
    t0 = time()
    Y = LocallyLinearEmbedding(n_neighbors, n_components,
                                        eigen_solver='auto',
                                        method=method).fit_transform(X)
    t1 = time()
    print("%s: %.2g sec" % (methods[i], t1 - t0))
    ax[1,i].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
    ax[1,i].set_title("%s (%.2g sec)" % (labels[i], t1 - t0))



#plt.figure()
#plt.scatter(np.arange(zmax.shape[0]),labelsm.T,c=colors[zmax], marker='o',s=50)
#plt.ylabel('Manual cell Label')
#plt.xlabel('Cell No')
#leg = ['GMM clusters']
#leg = leg[:len(set(zmax))]
#plt.legend(leg)