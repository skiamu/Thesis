#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 21:29:42 2018
script for retriving data from yahoo finance and get started with pandas
@author: andrea
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn import mixture
import scipy.io
from scipy.stats import norm

# import data from the .mat file
Returns = scipy.io.loadmat('AssetClassReturns.mat') # it's a dictionary
Returns = Returns['Returns']

# define an instance of class GaussianMixture
X = mixture.GaussianMixture(n_components=2,covariance_type='full',
                            n_init=50)
X.fit(Returns)
X.covariances_
X.weights_
X.means_

class GM(mixture.GaussianMixture):
    def __init__(self,X,**args):
        '''calls the super constructor and then calibrate the model '''
        mixture.GaussianMixture.__init__(self,**args)
        self.fit(X)
        self.M = X.shape[1] # asset allocation dimension
        
    def pf(self,z,u,x):
        '''method for computing the probability density function needed
        in the ODAA algorithm
        INPUT:
            z = points where the density is computed [vector]
            u = asset allocation vector [vector]
            x = realization portfolio value [scalar]
        OUTPUT:
            f = density computed in z'''
        f = 0
        for i in range(len(self.means_)):
            mu = x * (1+np.dot(u,self.means_[i]))
            sigma = np.sqrt(x**2 * self.covariances_[i].dot(u).dot(u))
            f = f + self.weights_[i] * norm.pdf(z,mu,sigma)
            return f
    def VaR_cons(self,VaR,alpha=0.01,u):
        W = 0
        for i in range(len(self.means_)):
            W = W + self.weights_[i] * self.covariances_[i]
            for j in range(i):
                W = W + self.weights_[i] * self.weights_[j] * \
                np.dot((self.means_[i]-self.means_[j]).reshape((self.M,1)),
                       self.means_[i]-self.means_[j])
        sigmaMax = VaR / norm.ppf(1-alpha)
        cons = {'type' : 'ineq',
                'fun' : sigmaMax-np.sqrt(W.dot(u).dot(u)),
                'jac' : W.dot(u)/np.sqrt(W.dot(u).dot(u))}
        return cons
    

Y = GM(X=Returns,n_components=2)
Y.means_


