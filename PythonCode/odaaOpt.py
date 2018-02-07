# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 16:19:02 2018

@author: andrea
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
from pyOpt import Optimization
from pyOpt import SLSQP
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
    for k in range(N-1,N-2,-1):
        print k # print current iteration
        u0 = np.array([1, 0, 0])
        dimXk = len(X[k])
        Uk = np.zeros((dimXk,model.M))
        Jk = np.zeros(dimXk)
        int_domain = np.arange(X[k+1][0],X[k+1][-1],eta)
        Jinterp = interp1d(X[k+1],J[k+1])(int_domain)
        for j in range(dimXk-1,dimXk-2,-1):
            #print j
            solveOpt(int_domain,Jinterp,X[k][j],model,u0)
            Uk[j] = res.x
            Jk[j] = res.fun
            u0 = Uk[j]
        U[k] = Uk
        J[k] = -Jk
        plt.figure
        plt.stackplot(X[k],Uk.T)
    return J,U
        
def solveOpt(int_domain,J,x,model,u0):
    def objfun(u,**kwargs):
        # 1) extract paraeters
        int_domain = kwargs['int_domain'] 
        J = kwargs['J'] 
        x = kwargs['x'] 
        model = kwargs['model'] 
        # 2) define objective function
        f = np.trapz(int_domain,J * model.pf(int_domain,u,x))
        g = [0]*2
        # 3) budget constraint 
        g[1] = u.sum() - 1
        # 4) VaR constarint
        W = model.W
        sigmaMax = model.VaR / norm.ppf(1-model.alpha)
        g[0] = -sigmaMax + np.sqrt(W.dot(u).dot(u))
        fail = 0
        return f,g,fail
    opt_prob = Optimization('test problem',objfun)
    opt_prob.addObj('f')
    opt_prob.addCon('budget const','e')    
    opt_prob.addCon('VaR const','i')
    opt_prob.addVarGroup('u',model.M,'c',lower=np.zeros(model.M),
                         upper=np.ones(model.M),value=u0)
    print opt_prob
    slsqp = SLSQP()
    slsqp.setOption('IPRINT',-1)
    slsqp(opt_prob,sens_type='FD',int_domain=int_domain,J=J,x=x,model=model)
    print opt_prob.solution(0)
    

      
    
    