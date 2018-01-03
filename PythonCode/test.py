#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 08:55:26 2017

@author: andrea
"""

import scipy.io # module for reading file .mat
import csv # for reading csv files
import numpy as np
import math
from scipy.optimize import minimize
import scipy.integrate as integrate
from scipy.stats import norm
import Calibration as cal
import matplotlib.pyplot as plt
# 1) read LIBOR file
v = list()
with open('LIBOR.csv','r') as csvfile:
    LIBOR = csv.reader(csvfile, delimiter=',')
    for row in LIBOR:
        v.append(row[1])
v = v[1:] 
LIBOR = [float(i) / 100  for i in v if i != '.']

# read stocks file
filemat = scipy.io.loadmat('RiskyAsset.mat')
RiskyAssetPrice = filemat['S']
RiskyAssetReturns = RiskyAssetPrice[1:] / RiskyAssetPrice[:-1] - 1


J = 0.07
L = cal.DiscretePrice(np.copy(RiskyAssetPrice),J)
D = L[0]
DiscreteS = L[1]
#xx = np.arange(len(D))
#plt.plot(xx,RiskyAssetPrice,'b',xx,DiscreteS,'r')

# optimization
bnds = ((0,1),(0,math.inf)) # upper and lower bounds
x0_p = np.random.uniform(0,1,100) # initial points
x0_lambda = 20 * np.random.uniform(0,1,100) # initial points
x = np.zeros(shape=(len(x0_p),2))
f = np.zeros(len(x0_p))

for i in range(len(x0_p)):
   res =  minimize(cal.objfun, [x0_p[i],x0_lambda[i]], args=((D,1/252,-1),),
                   method='SLSQP',bounds=bnds,options={'disp': False})
   f[i] = -res.fun
   x[i][:] = res.x
f_star = max(f)
x_star = x[np.argmax(f)]  
   


# integration test
def integrand(x, a, b):
    return a*x**2 + b
a = np.array([1,4])
a.shape = (2,1)
b = a 
x = np.linspace(0,1,20)
y = np.zeros(shape = (2,20))
y = integrand(x,a,b)

# integra lungo le righe
I = integrate.simps(y,x)
I
    

t = np.linspace(0,1,20)
