#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:52:04 2018

@author: andrea
"""
import numpy as np
import scipy.io
from returnsmodel import GM
import random
import odaaOpt2 

random.seed(1992)
# import data from the .mat file
Returns = scipy.io.loadmat('AssetClassReturns.mat') # it's a dictionary
Returns = Returns['Returns']
# define object GM
model = GM(X=Returns,n_components=2,n_init=500,tol=1e-5,
           reg_covar=0,warm_start=True,VaR=0.07/2,alpha=0.01)
# define ODAA parameters
N = 104;theta = 0.07; LB = 0.5; UB = 1.9
eta = 1e-3
X = [0]*(N+1)
for k in range(1,N):
    X[k] = np.arange(LB,UB,eta)
X[0] = 1; X[-1] = np.arange((1+theta)**2,UB,eta)    


J,U = odaaOpt2.ODAAalgorithmTD(N,X,model)
