# -*- coding: utf-8 -*-
"""
Created on Mon May 15 17:00:26 2017

@author: zubarei1
"""

import numpy as np
from time import time
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import entropy, kurtosis, skew
from scipy.special import cbrt
#from sklearn.linear_model import ARDRegression

from scipy.io import loadmat
#from bayespy.nodes import Poisson
path = '/m/home/home6/62/zubarei1/data/Desktop/elinas_APs/results/'
import os
os.chdir(path)
whiten=False
a = loadmat('AllCellsUpTo040517.mat')#loadmat('abfdata_27-Mar-2017_biocytin.mat')#
b = loadmat('BiocytinCells_04-May-2017.mat')
rf = []

sweep_data2 = a['X_2sw']

sweep_data = a['sweeps'][:,0]
spike_data = a['spikes']
cell_data = a['cells']
names = a['abflist']
#names = [np.array(n[0]) for n in names[0]]
X2sw = a['X_2sw']
X2swm = X2sw#18, 20, 21
x_ticks = ['Tau', 'ADP,mV','ADP,ms', 'Adapt,%','Step->NS1' ,'Step->NS2','Step->Delay1','Step->Delay2' ]

#GET TEST SET
spike_data_t = b['spikes']
sweep_data_t = b['sweeps'][:,0]
X2swt = b['X_2sw']#18, 20, 21
namest = b['abflist']
namest = [n[0] for n in namest[0]]



def poisson_regr(x,y,lamb=0):
    x = np.atleast_2d(x).T
    #x = np.vstack([np.ones(x.shape[-1]),x]).T
    y = np.log(np.atleast_2d(y)).T 
    XtX_lamb = x.T.dot(x) + lamb
    XtY = x.T.dot(y)
    b = np.squeeze(np.linalg.solve(XtX_lamb, XtY))
    return np.exp(b)
    
def lin_regr(x,y,lamb=0):
    x = np.atleast_2d(x)
    x = np.vstack([np.ones(x.shape[-1]),x]).T
    y = np.atleast_2d(y).T    
    XtX_lamb = x.T.dot(x) + lamb
    XtY = x.T.dot(y)
    b = np.squeeze(np.linalg.solve(XtX_lamb, XtY))
    return b
    
def gamma_regr(x,y):
    #x = np.vstack([np.ones(len(x)),np.log(np.array(x))]).T
    gamma_model = sm.GLM(y, x, family=sm.families.Gamma(link=sm.families.links.log))
    gamma_results = gamma_model.fit()
    bg = gamma_results.params
    return(bg)
 
ptr = []
ptst = []
rf = []   
"""Collect features"""    
for cellno, spk_sweeps in enumerate(spike_data):
    spk_sweeps = [spk for spk in spk_sweeps if np.any(spk)]
    predictor = sweep_data[cellno][:,1]
    ptr+=list(predictor)
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
    rf.append(regr_features)
    
rf = np.vstack(rf)
X = np.hstack([X2swm,rf])

rft = []   

"""Collect test set features"""    
for cellno, spk_sweeps in enumerate(spike_data_t):
    spk_sweeps = [spk for spk in spk_sweeps if np.any(spk)]
    predictor = sweep_data_t[cellno][:,1]
    ptst+=list(predictor)
    
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




"""Remove samples with bad/mising values"""
nso = np.array([np.shape(X)[0],np.shape(Xt)[0]])
nanrow, nancol = np.where(np.isnan(X))
if len(nanrow) !=0:
    print(np.unique(nanrow))
X = np.delete(X,np.unique(nanrow),0)
#filler = np.nanmean(X[:,nancol],0)[np.newaxis,...]
#X[[nanrow],[nancol]] = filler + 1e-3*np.random.randn(*filler.shape)
names = np.delete(names, np.unique(nanrow))
nanrow, nancol = np.where(np.isnan(Xt))
if len(nanrow) !=0:
    print(np.unique(nanrow))
Xt = np.delete(Xt,np.unique(nanrow),0)
#filler = np.nanmean(Xt[:,nancol],0)[np.newaxis,...]
#Xt[[nanrow],[nancol]] = filler + 1e-3*np.random.randn(*filler.shape)

namest = np.delete(namest, np.unique(nanrow))

"""Log-transform"""
logcols = np.array([1,2,5,6,8,17,18,19,20,21,22,24])#
X[:,logcols] = np.log(X[:,logcols])
Xt[:,logcols] = np.log(Xt[:,logcols])
cbrt_cols = np.array([9,10,12,23])
X[:,cbrt_cols] = cbrt(X[:,cbrt_cols])
Xt[:,cbrt_cols] = cbrt(Xt[:,cbrt_cols])


from sklearn.preprocessing import MinMaxScaler, RobustScaler
scaler =  RobustScaler()#MinMaxScaler()#
#X = scaler.fit_transform(X) 
#Xt = scaler.transform(Xt)
#
"""Drop Outlier Samples - POSITIVE"""
outcols = [2,21,22,23,25]
for outc in outcols:
    outrow = np.where(X[:,outc] > np.percentile(X[:,outc],95)*2.)[0]
    if len(outrow)!=0:
        print(outrow)
        X = np.delete(X,outrow,0)
        names = np.delete(names, np.unique(outrow))
    outrow = np.where(Xt[:,outc] > np.percentile(Xt[:,outc],95)*2.)[0]
    if len(outrow)!=0:
        print(outrow)
        Xt = np.delete(Xt,outrow,0)
        namest = np.delete(namest, np.unique(outrow))
#"""Drop Outlier Samples - NEGATIVE"""
#outcolsn = [24]
#for outc in outcolsn:
#    outrow = np.where(X[:,outc] < np.percentile(X[:,outc],5)*2)[0]
#    if len(outrow)!=0:
#        print(outrow)
#        X = np.delete(X,outrow,0)
#        names = np.delete(names, np.unique(outrow))
#    outrow = np.where(Xt[:,outc] < np.percentile(Xt[:,outc],5)*2)[0]
#    if len(outrow)!=0:
#        print(outrow)
#        Xt = np.delete(Xt,outrow,0)
#        namest = np.delete(namest, np.unique(outrow))

"""Drop Bad Features"""
old_ind = np.arange(X.shape[-1])
badfc = [0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,19,20,21,22,24,25,26,27]#21,22
X = np.delete(X,badfc,1)
Xt = np.delete(Xt,badfc,1)
old_ind = np.delete(old_ind,badfc)

"""Scale & Sort for RFE"""
X = scaler.fit_transform(X) 
Xt = scaler.transform(Xt)
#moments_train = np.vstack([X.mean(0),X.var(0), skew(X,0),kurtosis(X,0)])
#moments_test = np.vstack([Xt.mean(0),Xt.var(0), skew(Xt,0),kurtosis(Xt,0)])
#mcc = np.corrcoef(moments_train.T)
#fcc = np.corrcoef(X.T)
#forder = np.argsort(moments_train[-1,:])
#X = X[:,forder]
#Xt = Xt[:,forder]
#old_ind = old_ind[forder]

#
#
#for i,x in enumerate(X.T):
#    plt.hist(x,50); plt.hist(Xt[:,i],50,color='r',alpha=0.5)
#    plt.show()
#    print('TRAIN: %i, mean= %f3., var= %f3., skew= %f3., kurtosis= %f3., ' %(old_ind[i], np.mean(x), np.var(x), skew(x), kurtosis(x)))
#    print('TEST: %i, mean= %f3., var= %f3., skew= %f3., kurtosis= %f3., ' %(old_ind[i], np.mean(Xt[:,i]), np.var(Xt[:,i]), skew(Xt[:,i]), kurtosis(Xt[:,i])))
###    
    
nsn = np.array([np.shape(X)[0],np.shape(Xt)[0]])
print('Dropred TRAIN: %i, TEST: %i' %((nso-nsn)[0],(nso-nsn)[1]))




the_metric1 = []
the_metric2 = []
feature_inds = []
from sklearn import mixture
import itertools
##from scipy import linalg
##import matplotlib as mpl
##from sklearn.metrics import silhouette_score
lowest_bic_ever = np.infty
for j in np.arange(1,X.shape[-1]+1):
    Xrfe=X[:,:j]
    Xtrfe=Xt[:,:j]
    lowest_bic = np.infty
    lowest_bic_tt = np.infty
    bic = []
    lb = []
    silhuettes = []
    n_components_range = range(1, 7)
    cv_types = ['full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          random_state = 42,
                                          init_params='random',
                                          covariance_type=cv_type)
            gmm.fit(Xrfe)
            bic.append(gmm.bic(Xrfe)/Xrfe.shape[0])
            #labels = gmm.predict(X)
            ll = gmm.bic(Xtrfe)/Xtrfe.shape[0]
            lb.append(ll)
            #if n_components > 1:
            #silhuettes.append(silhouette_score(X,labels,metric='mahalanobis'))
            
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm_bic = gmm
            if lb[-1]+bic[-1] < lowest_bic_tt:
                lowest_bic_tt = lb[-1]+bic[-1]
                if lowest_bic_tt< lowest_bic_ever:
                    best_gmm = gmm
    the_metric1.append(lowest_bic/j)
    the_metric2.append(lowest_bic_tt/j)
#        
#   
#    bic = np.array(bic)
#    lb = np.array(lb)
#    #silhuettes = np.array(silhuettes)
#    color_iter = itertools.cycle(['turquoise'])
#    clf = best_gmm
#    bars = []
#    
#    # Plot the BIC scores
#    spl = plt.subplot(2, 1, 1)
#    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
#        xpos = np.array(n_components_range) + .2 * (i - 2)
#        bars.append(plt.bar(xpos, bic[i * len(n_components_range):
#                                      (i + 1) * len(n_components_range)],
#                            width=.2, color=color))
#    plt.xticks(n_components_range)
#    plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
#    plt.title('TRAIN BIC')
#    xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
#        .2 * np.floor(bic.argmin() / len(n_components_range))
#    plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
#    spl.set_xlabel('Number of components')
#    spl.legend([b[0] for b in bars], cv_types)
#    
#    # Plot the winner
#    splot = plt.subplot(2, 1, 2)
#    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
#        xpos = np.array(n_components_range) + .2 * (i - 2)
#        bars.append(plt.bar(xpos, lb[i * len(n_components_range):
#                                      (i + 1) * len(n_components_range)],
#                            width=.2, color=color))
#    plt.xticks(n_components_range)
#    plt.ylim([lb.min() * 1.01, lb.max()*1.01])
#    plt.title('TEST BIC')
#    xpos = np.mod(lb.argmin(), len(n_components_range)) + .65 +\
#        .2 * np.floor(lb.argmin() / len(n_components_range))
#    plt.text(xpos, lb.min()* 0.97 + .03 * lb.max(), '*', fontsize=14)
#    splot.set_xlabel('Number of components')
#    splot.legend([b[0] for b in bars], cv_types)
#    print(old_ind[j-1])
#    plt.show()
#
plt.plot(the_metric1)
plt.hold(True)
plt.plot(the_metric2,'r--')
print(best_gmm.weights_)
plt.show()

fin=True
if fin:
    ind = np.argsort(the_metric2)[0]
    Xf = X[:,:ind+1]
    Xft = Xt[:,:ind+1]
    lowest_bic = np.infty
    lowest_bic_tt = np.infty
    bic = []
    lb = []
    silhuettes = []
    n_components_range = range(1, 7)
    cv_types = ['diag','full']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a Gaussian mixture with EM
            gmm = mixture.GaussianMixture(n_components=n_components,
                                          random_state = 42,
                                          init_params='random',
                                          covariance_type=cv_type)
            gmm.fit(Xf)
            bic.append(gmm.bic(Xf)/Xf.shape[0])
            #labels = gmm.predict(X)
            ll = gmm.bic(Xft)/Xft.shape[0]
            lb.append(ll)
            #if n_components > 1:
            #silhuettes.append(silhouette_score(X,labels,metric='mahalanobis'))
            
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm_bic = gmm
            if lb[-1]+bic[-1] < lowest_bic_tt:
                lowest_bic_tt = lb[-1]+bic[-1]
                if lowest_bic_tt< lowest_bic_ever:
                    best_gmm = gmm
    bic = np.array(bic)
    lb = np.array(lb)
    color_iter = itertools.cycle(['turquoise', 'cornflowerblue'])
    bars = []
    # Plot the BIC scores
    spl = plt.subplot(2, 1, 1)
    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        xpos = np.array(n_components_range) + .2 * (i - 2)
        bars.append(plt.bar(xpos, bic[i * len(n_components_range):
                                      (i + 1) * len(n_components_range)],
                            width=.2, color=color))
    plt.xticks(n_components_range)
    plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
    plt.title('TRAIN BIC')
    xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
        .2 * np.floor(bic.argmin() / len(n_components_range))
    plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
    spl.set_xlabel('Number of components')
    spl.legend([b[0] for b in bars], cv_types)
    
    # Plot the winner
    splot = plt.subplot(2, 1, 2)
    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        xpos = np.array(n_components_range) + .2 * (i - 2)
        bars.append(plt.bar(xpos, lb[i * len(n_components_range):
                                      (i + 1) * len(n_components_range)],
                            width=.2, color=color))
    plt.xticks(n_components_range)
    plt.ylim([lb.min() * 1.01, lb.max()*1.01])
    plt.title('TEST BIC')
    xpos = np.mod(lb.argmin(), len(n_components_range)) + .65 +\
        .2 * np.floor(lb.argmin() / len(n_components_range))
    plt.text(xpos, lb.min()* 0.97 + .03 * lb.max(), '*', fontsize=14)
    splot.set_xlabel('Number of components')
    splot.legend([b[0] for b in bars], cv_types)
    print(old_ind[:ind+1])
    print(best_gmm.weights_)
    plt.show()

#
#
#
#
#
#
#
#






