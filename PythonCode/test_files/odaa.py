#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 17:11:14 2018

@author: andrea
"""

from scipy.optimize import minimize
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt


def ODAAalgorithmTD(N,X,model,VaR,alpha=0.01):
    ''' 
    INPUT:
        N = number of time steps
        X = discretized target sets [tuple dimension N+1]
        model = object containing model information and parameters
        VaR = 
        alpha =
    '''
    U = [0]*N 
    J = [0]*(N+1)
    J[-1] = np.ones(len(X[-1]))
    eta = 1e-4/2 # discretization step
    budget_cons = {'type':'eq','fun':lambda u : np.sum(u)-1,
                   'jac':lambda u : np.ones(model.M)}
    for k in range(N-1,N-2,-1):
        print k # print current iteration
        u0 = np.array([1, 0, 0])
        dimXk = len(X[k])
        Uk = np.zeros((dimXk,model.M))
        Jk = np.zeros(dimXk)
        int_domain = np.arange(X[k+1][0],X[k+1][-1],eta)
        Jinterp = interp1d(X[k+1],J[k+1])(int_domain)
        Bounds = ((0,1),)*model.M
        for j in range(dimXk-1,-1,-1):
            #print j
            res = minimize(obj_fun,u0,args=(X[k][j],int_domain,Jinterp,model,-1),
                           method = 'SLSQP', bounds = Bounds,
                           constraints = (model.VaR_cons(VaR,alpha),budget_cons),
                           options={'disp': False} )
            Uk[j] = res.x
            Jk[j] = res.fun
            u0 = Uk[j]
        U[k] = Uk
        J[k] = -Jk
        plt.figure
        plt.stackplot(X[k],Uk.T)
        #input('press enter to continue')
    return J,U
        
        
def obj_fun(u,x,int_domain,J,model,sign):
    f = np.trapz(x = int_domain, y = J * model.pf(int_domain,u,x))
    return sign * f
        
        
    
    
    