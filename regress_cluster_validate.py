# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:29:26 2016

@author: zubarei1
"""

import numpy as np
from time import time
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
#from sklearn.linear_model import ARDRegression

from scipy.io import loadmat
#from bayespy.nodes import Poisson
path = '/m/home/home6/62/zubarei1/data/Desktop/elinas_APs/results/'
import os
os.chdir(path)
whiten=False
a = loadmat('abfdata_27-Mar-2017_biocytin.mat')#loadmat('abfdata_17-Mar-2017.mat')#
rf = []

sweep_data2 = a['X_2sw']

sweep_data = a['sweeps'][:,0]
spike_data = a['spikes']
cell_data = a['cells']
names = a['abflist']
#names = [np.array(n[0]) for n in names[0]]
X2sw = a['X_2sw']
X2swm = X2sw[:,[5,9,10,23]]#18, 20, 21
x_ticks = ['Tau', 'ADP,mV','ADP,ms', 'Adapt,%','Step->NS1' ,'Step->NS2','Step->Delay1','Step->Delay2' ]

#GET TEST SET
b = loadmat('abfdata_27-Mar-2017_together.mat')
spike_data_t = b['spikes']
sweep_data_t = b['sweeps'][:,0]
X2swt = b['X_2sw'][:,[5,9,10,23]]#18, 20, 21
namest = b['abflist']
namest = [n[0] for n in namest[0]]



def poisson_regr(y,x,lamb=0):
    x = np.atleast_2d(x)
    x = np.vstack([np.ones(x.shape[-1]),x]).T
    y = np.log(np.atleast_2d(y)).T 
    XtX_lamb = x.T.dot(x) + lamb
    XtY = x.T.dot(y)
    b = np.squeeze(np.linalg.solve(XtX_lamb, XtY))
#    f = plt.figure()
#    model = np.dot(b,x.T)
#    plt.hold(True)
#    plt.semilogy(x[:,-1],model);plt.semilogy(x[:,-1],y,'go')
#    plt.show()
    return b
    

#def poisson_regr(x,y,lamb=0):
#    poisson_model = sm.Poisson(endog=y,exog=x)
#    pr = poisson_model.fit()
#    b = pr.params
#    model = np.exp(b*x)
#    plt.hold(True)
#    plt.plot(x,model);plt.plot(x,y,'go')
#    plt.show()
#    return b
    
def lin_regr(x,y,lamb=0):
    x = np.atleast_2d(x)
    x = np.vstack([np.ones(x.shape[-1]),x]).T
    y = np.atleast_2d(y).T    
    XtX_lamb = x.T.dot(x) + lamb
    XtY = x.T.dot(y)
    b = np.squeeze(np.linalg.solve(XtX_lamb, XtY))
#    model = np.dot(b,x.T)
#    plt.plot(model);plt.plot(y)
#    plt.show()
    return b
    
def gamma_regr(x,y):
    
    x = np.vstack([np.ones(len(x)),np.log(np.array(x))]).T
    # Instantiate a gamma family model with the default link function.
    gamma_model = sm.GLM(y, x, family=sm.families.Gamma(link=sm.families.links.log))
    gamma_results = gamma_model.fit()
    bg = gamma_results.params
#    model = np.dot(bg,x.T)
#    f = plt.figure()
#    plt.plot(np.exp(x[:,-1]),np.exp(model));plt.hold(True),plt.plot(np.exp(x[:,-1]),y,'ko')
#    plt.show()
#    print(bg)
    return(bg)
 

#def broken_stick(x,y,define=False):
#    #x = predictor
#    #y = count
#    """Define if the stick is broken"""
#    if define:
#        print('NOT IMPLENTED')
#    else:    
#        bs_ind = np.where(y==np.max(y))[0][-1]
#    x1 = x[:bs_ind+1]
#    y1 = y[:bs_ind+1]
#    x2 = x[bs_ind:]
#    y2 = y[bs_ind:]
#    
#    if len(y2) > 5:
#        broken= 1
#        b1 = poisson_regr(x1,y1)
#        b2 = poisson_regr(x2,y2)
#    else:
#        broken = 0
#        bs_ind = 0
#        b1 = poisson_regr(x,y)
#        b2 = np.zeros(2,)
#    model1 = np.dot(b1,np.vstack([np.ones(len(x1)),x1]))
#    model2 = np.dot(b2,np.vstack([np.ones(len(x2)),x2]))
#    if broken==0:
#        model1 = np.exp(model1)
#        model2 = np.exp(model2)
   # plt.plot(x1,np.exp(model1),'b');
   # plt.plot(x,(y),'k');
    #plt.plot(x2,np.exp(model2),'b')
    #plt.show()
    #print([b1[0]*np.exp(b1[1]),b2[0]*np.exp(b2[1])])
    #print([b1[1], b2[1]])
    #return np.hstack([broken,x[bs_ind],b1,b2])
    
    
   
rf = []   
"""Collect features"""    
for cellno, spk_sweeps in enumerate(spike_data):
    spk_sweeps = [spk for spk in spk_sweeps if np.any(spk)]
    predictor = sweep_data[cellno][:,1]
    
    """Gamma-Poisson distributed features"""
    first_lat= []
    count = []
  
    for sweep in spk_sweeps:
        if np.any(sweep):
            first_lat.append(sweep[0,0])
            count.append(sweep.shape[0])
    inp_nspk = gamma_regr(predictor,count)
    inp_delay = gamma_regr(predictor,first_lat)
    regr_features = np.hstack([inp_nspk,inp_delay])
    rf.append(regr_features)
    
rf = np.vstack(rf)
X = np.hstack([X2swm,rf])

rft = []   

"""Collect test set features"""    
for cellno, spk_sweeps in enumerate(spike_data_t):
    spk_sweeps = [spk for spk in spk_sweeps if np.any(spk)]
    predictor = sweep_data_t[cellno][:,1]
    
    """Gamma-Poisson distributed features"""
    first_lat= []
    count = []
  
    for sweep in spk_sweeps:
        if np.any(sweep):
            first_lat.append(sweep[0,0])
            count.append(sweep.shape[0])
    inp_nspk = poisson_regr(predictor,count)
    inp_delay = gamma_regr(predictor,first_lat)
    regr_features = np.hstack([inp_nspk,inp_delay])
    rft.append(regr_features)
    
rft = np.vstack(rft)
Xt = np.hstack([X2swt,rft])




"""Remove samples with bad/miisng values"""
nanrow, nancol = np.where(np.isnan(X))
X = np.delete(X,np.unique(nanrow),0)
names = np.delete(names, np.unique(nanrow))
outrow = np.where(X[:,3] > np.percentile(X[:,3],95)*2)
X = np.delete(X,outrow,0)
names = np.delete(names, np.unique(outrow))
#
#"""Remove samples with bad/miisng values"""
nanrow, nancol = np.where(np.isnan(Xt))
Xt = np.delete(Xt,np.unique(nanrow),0)
namest = np.delete(namest, np.unique(nanrow))
outrow = np.where(Xt[:,3] > np.percentile(X[:,3],95)*2)
Xt = np.delete(Xt,outrow,0)
namest = np.delete(namest, np.unique(outrow))
#
#
#"""Scale/Whiten"""
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
scaler = MinMaxScaler()
X = scaler.fit_transform(X) 
Xt = scaler.transform(Xt)

if whiten:
    from sklearn.decomposition import PCA
    pca = PCA(n_components=0.99, whiten=True)
    X = pca.fit_transform(X)
    Xt = pca.transform(Xt)

#from sklearn.manifold import TSNE, MDS, SpectralEmbedding, LocallyLinearEmbedding
from bayespy.nodes import Dirichlet, Categorical
from bayespy.nodes import  GaussianARD, Gamma, Gaussian, Wishart
from bayespy.nodes import Mixture
from bayespy.inference import VB

"""Train GMM MODEL"""
plt.close('all')
sc_alpha = 1e-7
sc_m0 = 1
sc_gam_rate = 1e-2
sc_gam_shape = 1e-2
n, d = X.shape
k = 8
a0 = sc_alpha*np.ones(k)/k
alpha = Dirichlet(a0, name='alpha')

#WITH ARD
Z = Categorical(alpha, plates=(n,1), name='Z')
m0 = GaussianARD(0.5,sc_m0,shape=(d,), plates=(k,), name='m0')
lam = Gamma(sc_gam_shape,sc_gam_rate,plates=(k,d,), name='lam')
Y = Mixture(Z, GaussianARD, m0, lam, name='Y',cluster_plate=-2)

#WITHOUT ARD
#Z = Categorical(alpha, plates=(n,), name='Z')
#m0 = Gaussian(np.zeros(d), sc_m0*np.identity(d), plates=(k,), name='m0')
#lam = Wishart(d, sc_gam_rate*np.identity(d), plates=(k,), name='lam')
#Y = Mixture(Z, Gaussian, m0, lam, name='Y')

Z.initialize_from_random()
Q = VB(Y, m0, lam, Z, alpha)
Y.observe(X)
Q.update(repeat=1000)
means0 = Q['m0'].get_moments()[0]
precs0 = Q['lam'].get_moments()[0]
alphas = Q['alpha'].get_moments()[0]
Zs = Q['Z'].get_moments()[0]
ass0 = np.exp(alphas)
print(ass0)





"""VISUALIZE"""
ass1=[]
means = []
precs = []
Zm = []
ignoresmall = 0.025
Zc = np.squeeze(Zs)
for i, a0 in enumerate(np.exp(alphas)):
    if a0 >ignoresmall:
        means.append(means0[i,...])
        precs.append(precs0[i,...])
        ass1.append(ass0[i])
        Zm.append(Zc[:,i])
ind = np.argsort(np.array(ass1))[::-1]
#act_comp = np.where(np.exp(alphas)>1e-6)[0]

#Zc = Zc[...,act_comp]
Zc = np.vstack(Zm).T[...,ind]
zmax = np.argmax(Zc,-1)
ass1 = np.array(ass1)[ind]
means = np.array(means)[ind]
precs = np.array(precs)[ind]
leg = ['Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Ccluster8']
colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
leg = leg[:len(ass1)]

plt.figure()
plt.hold(True)
[plt.errorbar(np.arange(X.shape[-1]),means[i],yerr=np.sqrt(1./precs[i]),marker='s', linestyle='None', color=colors[i],markersize=np.sqrt(ass1[i])*20) for i in range(len(means)) if ass1[i]>ignoresmall]
plt.legend(leg)#plt.errorbar(np.arange(7),means[1],yerr=1./precs[1],marker='o',color='r');
plt.xlim(-1,len(means[0]))
plt.xticks(np.arange(X.shape[-1]),x_ticks,fontsize=14)
#plt.xlabel(*x_ticks)

plt.show()

#np.savez('1703_data_biocyt', zmax=zmax, Z=Zc, proportion=ass1, means=means,stds=np.sqrt(1./precs))
np.savetxt('2703_data_biocytin.csv', np.hstack([names[...,None],zmax[...,None],Zc]),fmt='%s,%f,%f,%f')
#
#from sklearn.metrics import silhouette_samples, silhouette_score, calinski_harabaz_score
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import numpy as np
#from sklearn.cluster import KMeans, AgglomerativeClustering, MeanShift, AffinityPropagation
#silhouette_avg = silhouette_score(X, zmax)
#print(silhouette_avg)
#chs = calinski_harabaz_score(X, zmax)
#print(chs)  
#
#f,ax = plt.subplots(2,3)
#n_neighbors = 5
#n_components = 3
#
#
#colors = np.array(['r', 'g', 'b', 'y', 'c', 'm', 'k', 'pink'])#, [0.5,0.5,0.5], [0.75, 0.75,0]])
#t0 = time()
#mds = MDS(n_components, max_iter=1000, n_init=1)
#Y = mds.fit_transform(X)
#t1 = time()
#print("MDS: %.2g sec" % (t1 - t0))
##ax = fig.add_subplot(258)
##a[].figure()
#ax[0,0].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
#ax[0,0].set_title("MDS (%.2g sec)" % (t1 - t0))
#
#t0 = time()
#se = SpectralEmbedding(n_components=n_components, n_neighbors=n_neighbors)
#Y = se.fit_transform(X)
#t1 = time()
#print("SpectralEmbedding: %.2g sec" % (t1 - t0))
#ax[0,1].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
#ax[0,1].set_title("SpectralEmbedding (%.2g sec)" % (t1 - t0))
#
#Xp = X - X.min(0)
#t0 = time()
#tsne = TSNE(n_components=n_components, init='pca', random_state=42)
#Y = tsne.fit_transform(Xp)
#t1 = time()
#print("t-SNE: %.2g sec" % (t1 - t0))
#ax[0,2].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
#ax[0,2].set_title("t-SNE (%.2g sec)" % (t1 - t0))
#
#methods = ['standard', 'modified', 'ltsa']
#labels = ['LLE', 'Modified LLE', 'LTSA']
#
#for i, method in enumerate(methods):
#    t0 = time()
#    Y = LocallyLinearEmbedding(n_neighbors, n_components,
#                                        eigen_solver='auto',
#                                        method=method).fit_transform(X)
#    t1 = time()
#    print("%s: %.2g sec" % (methods[i], t1 - t0))
#    ax[1,i].scatter(Y[:, 0], Y[:, 1], c=colors[zmax], cmap=plt.cm.Spectral)
#    ax[1,i].set_title("%s (%.2g sec)" % (labels[i], t1 - t0))
##
##
#
#
#
##
#
#
#n_clusters = 3
#range_n_methods = [MeanShift(),
##                    AffinityPropagation(damping=0.5),
##                    AffinityPropagation(damping=0.75),
##                    AffinityPropagation(damping=0.95),
##                   AgglomerativeClustering(n_clusters=n_clusters,
##                   affinity='cosine', linkage='average'),
#                    AgglomerativeClustering(n_clusters=n_clusters,
#                    affinity='manhattan', linkage='complete'),
#                   AgglomerativeClustering(n_clusters=n_clusters,
#                   affinity='cosine', linkage='complete')
#                   ]
#                   
#                   #, random_state=42)]
#
#for clusterer in range_n_methods:
#    # Create a subplot with 1 row and 2 columns
#    #fig, (ax1, ax2) = plt.subplots(1, 2)
#    fig = plt.figure()
#    fig.set_size_inches(6, 6)
#
#    # The 1st subplot is the silhouette plot
#    # The silhouette coefficient can range from -1, 1 but in this example all
#    # lie within [-0.1, 1]
#    plt.xlim([-0.1, 1])
#    # The (n_clusters+1)*10 is for inserting blank space between silhouette
#    # plots of individual clusters, to demarcate them clearly.
#    plt.ylim([0, len(X) + (n_clusters + 1) * 10])
#
#    # Initialize the clusterer with n_clusters value and a random generator
#    # seed of 10 for reproducibility.
#    #clusterer = KMeans(n_clusters=n_clusters, random_state=10)
#    #clusterer = MeanShift()
#    #clusterer = AffinityPropagation(damping=0.99)
#    #AgglomerativeClustering(n_clusters=n_clusters,
#                #                        affinity='cosine', linkage='average')#, random_state=42)
#    cluster_labels = clusterer.fit_predict(X)
#
#    # The silhouette_score gives the average value for all the samples.
#    # This gives a perspective into the density and separation of the formed
#    # clusters
#    silhouette_avg = silhouette_score(X, cluster_labels)
#    print("For n_clusters =", clusterer,
#          "The average silhouette_score is :", silhouette_avg)
#
#    # Compute the silhouette scores for each sample
#    sample_silhouette_values = silhouette_samples(X, cluster_labels)
#
#    y_lower = 10
#    for i in range(n_clusters):
#        # Aggregate the silhouette scores for samples belonging to
#        # cluster i, and sort them
#        ith_cluster_silhouette_values = \
#            sample_silhouette_values[cluster_labels == i]
#
#        ith_cluster_silhouette_values.sort()
#
#        size_cluster_i = ith_cluster_silhouette_values.shape[0]
#        y_upper = y_lower + size_cluster_i
#
#        color = cm.spectral(float(i) / n_clusters)
#        plt.fill_betweenx(np.arange(y_lower, y_upper),
#                          0, ith_cluster_silhouette_values,
#                          facecolor=color, edgecolor=color, alpha=0.7)
#
#        # Label the silhouette plots with their cluster numbers at the middle
#        plt.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
#
#        # Compute the new y_lower for next plot
#        y_lower = y_upper + 10  # 10 for the 0 samples
#
#    plt.title("The silhouette plot for the various clusters.")
#    plt.xlabel("The silhouette coefficient values")
#    plt.ylabel("Cluster label")
#
#    # The vertical line for average silhouette score of all the values
#    plt.axvline(x=silhouette_avg, color="red", linestyle="--")
#
#    plt.yticks([])  # Clear the yaxis labels / ticks
#    plt.xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
#
#    # 2nd Plot showing the actual clusters formed
#    colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
#    plt.show()
##    
