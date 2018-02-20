# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:29:01 2018

@author: andrea
"""
from __future__ import division # import division operator from 3.X
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from scipy.integrate import trapz,simps
from scipy.interpolate import interp1d

# =========================================================================
#   CLASS DEFINITION
# ========================================================================
class BasicModel:
    def __init__(self,J_jump,r,dt,S):
        self.J_jump = J_jump
        self.r = r
        self.dt = dt
        self.p,self.Lambda = self.Calibration(S)
    
    def Calibration(self,S):
        D,discreteS = DiscretePrice(S,self.J_jump)
        alpha = len(D[D==0])
        beta = len(D[D==1])
        gamma = len(D[D==-1])
        p = beta / (beta + gamma)
        Lambda = -np.log(alpha / len(D)) / self.dt
        return p,Lambda
    
    def pf(self,z,u,x):
        f = np.zeros(z.shape)
        csi = x * self.J_jump * u
        idx1 = z >= x + csi
        idx2 = z >= x - csi
        f[idx1] = self.p*((z[idx1]-csi)/x)**(-(self.Lambda+self.r)/self.r)
        f[idx2] = f[idx2]+(1-self.p)*((z[idx2]+csi)/x)**(-(self.Lambda+self.r)/self.r)
        return self.Lambda/(self.r*x)*f
        
    def Simulate_rv(self,Nsim,Nstep):
        tau = np.random.exponential(1/self.Lambda,(Nsim,Nstep))
        Binomial = -np.ones((Nsim,Nstep))
        Binomial[np.random.uniform(0,1,(Nsim,Nstep))<self.p] = 1
        return tau, Binomial

class GBMmodel:
    def __init__(self,J_jump,r,dt,LogReturns):
        self.J_jump = J_jump
        self.r = r
        self.dt = dt
        self.mu,self.sigma = self.Calibration(LogReturns)
        mu_tilde = self.mu - 0.5*self.sigma**2
        self.p = (np.exp(2*mu_tilde*self.J_jump/self.sigma**2) - 1) / \
        (2*np.sinh(2*mu_tilde*self.J_jump/self.sigma**2))
        
    def Calibration(self,LogReturns):
        sample_var = np.var(LogReturns)
        sample_mean = np.mean(LogReturns)
        sigma = np.sqrt(sample_var/self.dt)
        mu = sample_mean / self.dt + 0.5*sigma**2
        return mu,sigma
     
    def pf(self,z,u,x):
        f = np.zeros(z.size)
        csi = x*self.J_jump*u
        idx1 = z > x + csi + np.finfo(float).eps
        idx2 = z > x - csi + np.finfo(float).eps
        mu_tilde = self.mu - 0.5*self.sigma**2
        if z[idx1].size is not 0:
            f[idx1] = self.p * Gamma(self,(z[idx1]-csi)/x)
        if z[idx2].size is not 0:
            f[idx2] = f[idx2] + (1-self.p) * Gamma(self,(z[idx2]+csi)/x)
        f = 2*np.cosh(mu_tilde*self.J_jump/self.sigma**2)/(self.r*x)*f
        return f

    def Simulate_rv(self,Nsim,Nstep):
        mu_tilde = self.mu-.5*self.sigma**2
        F_tau = lambda t : 1-2*np.cosh(mu_tilde*self.J_jump/self.sigma**2)*\
            Kappa(self,t)
        t = np.arange(1e-3,4,1e-5)
        u,idxUnique = np.unique(F_tau(t),return_index=True)
        u_tilde = np.random.uniform(0,1,(Nsim,Nstep))
        tau = interp1d(u,t[idxUnique])(u_tilde)
        Binomial = -np.ones((Nsim,Nstep))
        Binomial[np.random.uniform(0,1,(Nsim,Nstep))<self.p] = 1
        return tau, Binomial
        
        
# ========================================================================
#   AUXILIARY FUNCTIONs
# ========================================================================       
        
def DiscretePrice(S,J):
    """
    INPUT:
        S = prive vector
        J = jump size
    OUTPUT:
        D = 
        DiscreteS = 
    """
    N = len(S)
    DiscreteS = S # initialize discrete path vector
    D = np.zeros(N) # initialize jump sign vector
    Baseline = S[0]
    for i in range(1,N):
        r = S[i]/Baseline-1
        if np.abs(r) > J:
            Baseline = Baseline*(1+np.sign(r)*J) # update baseline
            DiscreteS[i] = Baseline # update discrete path
            D[i] = np.sign(r) # update jump sign vector
        else:
            DiscreteS[i] = DiscreteS[i-1]
            D[i] = 0
    return D,DiscreteS   

def Gamma(model,y):
    eps = 1e-8
    mu_tilde = model.mu - 0.5*model.sigma**2
    t_min = 1e-4 # minumum value such that there's no memory problem
    t_vec = np.log(y)/model.r
    t = np.min(t_vec)
    n = len(y) # save the length of the original y
    Flag = False
    if t <= t_min and n > 1:
        idx = t_vec > t_min # indexes of elements greater than t_min
        t = np.min(t_vec[idx]) # compute new t
        y = y[idx] # compute new y without the small elements
        Flag = True
    Ntrunc = np.sqrt(np.max([1,-8*model.J_jump**2/(np.pi**2*model.sigma**2*t) \
        *(np.log((np.pi**3*model.sigma**2*t*eps)/(16*model.J_jump**2))- \
        model.J_jump*mu_tilde/model.sigma**2)]))
    N = np.arange(1,Ntrunc+1)
    N = N.reshape((len(N),1)) # column vector
    f = model.sigma**2*np.pi/(4*model.J_jump**2)*(np.sum(N*(-1)**(N+1)*y**\
        (-(mu_tilde**2/(2*model.sigma**2)+(model.sigma**2*N**2*np.pi**2)/\
        (8*model.J_jump**2))/model.r - 1)*np.sin(np.pi*N/2),axis=0))
    #f[f<0] = 0    
    #assert all(f>0),'there are negative elements' 
    if Flag:
        g = np.zeros(n) # set to zero elements associated with small t values
        g[idx] = f
        return g
    return f    
    
def Kappa(model,y):
    eps = 1e-8
    mu_tilde = model.mu - 0.5*model.sigma**2
    t = np.min(y)
    Ntrunc = np.sqrt(np.max([1,-8*model.J_jump**2/(np.pi**2*model.sigma**2*t) \
        *(np.log((np.pi**3*model.sigma**2*t*eps)/(16*model.J_jump**2))- \
        model.J_jump*mu_tilde/model.sigma**2)]))
    N = np.arange(1,Ntrunc+1)
    N = N.reshape((len(N),1)) # column vector
    f = model.sigma**2*np.pi/(4*model.J_jump**2)*(np.sum(N*(-1)**(N+1)/\
        (mu_tilde**2/(2*model.sigma**2)+(model.sigma**2*N**2*np.pi**2)/\
        (8*model.J_jump**2))*\
        np.exp(-(mu_tilde**2/(2*model.sigma**2)+(model.sigma**2*N**2*np.pi**2)/\
        (8*model.J_jump**2))*y)*np.sin(np.pi*N/2),axis=0))
    return f    
    

#=========================================================================
#   TEST SCRIPT
#=========================================================================    

if __name__ == '__main__':
    modelUsed = 'basic'
    Data = scipy.io.loadmat('Data.mat') # it's a dictionary
    Returns = Data['Returns']
    S = Data['S']
    J_jump = 0.07
    r = 0.055
    dt = 1/252
    if modelUsed is 'basic':
        model = BasicModel(J_jump,r,dt,S)
        print model.p,model.Lambda,dt
        xx = np.arange(.5,1.9,1e-6)
        u = 0.0;x = 1
        print trapz(x = xx,y = model.pf(xx,u,x))
        plt.plot(xx,model.pf(xx,u,x),'.-')
        a,b = model.Simulate_rv(10,1)
    elif modelUsed is 'ext1':
        LogReturns = np.log(1+Returns)
        model = GBMmodel(J_jump,r,dt,LogReturns)
        #print model.p,model.mu,model.sigma
        xx = np.arange(.5,1.9,1e-4/3)
        u = 0.8;x = 1
        print trapz(x = xx,y = model.pf(xx,u,x))
        plt.plot(xx,model.pf(xx,u,x),'.-')
        a,b = model.Simulate_rv(10,1)
    
    
