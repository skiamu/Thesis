# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:43:36 2018

@author: andrea
"""
from __future__ import division # import division operator from 3.X
import scipy.io
from models import BasicModel,GBMmodel
import odaaOptED
import numpy as np
import shelve
import time
import matplotlib.pyplot as plt

# ========================================================================
#   RUN SCRIPT
# ========================================================================
modelUsed = 'ext1'
Data = scipy.io.loadmat('Data.mat') # it's a dictionary
Returns = Data['Returns']
S = Data['S']
J_jump = 0.03
r = 0.055
dt = 1/252
if modelUsed is 'basic':
    model = BasicModel(J_jump,r,dt,S)
elif modelUsed is 'ext1':
    LogReturns = np.log(1+Returns)
    model = GBMmodel(J_jump,r,dt,LogReturns)
N = 10;theta = 0.07; LB = 0.5; UB = 1.9
eta = 1e-3/2
X = [0]*(N+1)
for k in range(1,N):
    X[k] = np.arange(LB,UB+eta,eta)
X[0] = np.array([1]); X[-1] = np.arange((1+theta)**2,UB,eta)    

# run ODAA algorithm
start = time.clock()
J,U = odaaOptED.ODAAalgorithmED(N,X,model)
t = time.clock()-start
# saving the results is a shelf
db = shelve.open('ext1') # open a shelve
db['J'] = J 
db['U'] = U
db.close() # close the shelve

i = 1
plt.figure()
for k in range(N-1,0,-2):
    plt.subplot(3,2,i)
    plt.stackplot(X[k],U[k])
    plt.title('k = ' + str(k))
    i+=1

















