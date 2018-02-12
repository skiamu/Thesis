#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 17:07:21 2018

@author: andrea

This modele contains classes used to described asset class returns.
The models considered are:
    1) Gaussian (G)
    2) Gaussian Mixture (GM)
    3) Generalized Hyperbolic (GH)
"""
import numpy as np
from scipy.stats import norm
from sklearn import mixture

class GM(mixture.GaussianMixture):
    def __init__(self,X,alpha,VaR,**args):
        '''calls the super constructor and then calibrate the model '''
        mixture.GaussianMixture.__init__(self,**args)
        self.fit(X)
        self.M = X.shape[1] # asset allocation dimension
        self.alpha = alpha 
        self.VaR = VaR
        self.W = self.compute_W()
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
            sigma = np.sqrt(x**2 * (self.covariances_[i].dot(u).dot(u)))
            f = f + self.weights_[i] * norm.pdf(z,mu,sigma)
        return f
    
    def VaR_cons(self,VaR,alpha=0.01): # needed only when scipy.optimize is used
        W = self.compute_W()
        sigmaMax = VaR / norm.ppf(1-alpha)
        cons = {'type' : 'ineq',
                'fun' : lambda u : sigmaMax - np.sqrt(W.dot(u).dot(u)),
                'jac' : lambda u : (W.dot(u)) / np.sqrt(W.dot(u).dot(u))}
        return cons
        
    def compute_W(self):
        W = 0
        for i in range(len(self.means_)):
            W = W + self.weights_[i] * self.covariances_[i]
            for j in range(i):
                W = W + self.weights_[i] * self.weights_[j] * \
               (self.means_[i]-self.means_[j]).reshape((self.M,1)) * \
                       (self.means_[i]-self.means_[j])
        return W



