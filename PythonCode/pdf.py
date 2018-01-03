#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 09:54:00 2017

@author: andrea
"""
import scipy.integrate as integrate
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import time
# define parameters
a = 0.1670; b = 0.0280; sigma = 0.0160; p = 0.5745; lambda_p = 5.7820
r0 = 0.0550; J_jump = 0.07

# define auxiliary function
def mu_tilde(t):
    return ((r0-b)*(1-np.exp(-a*t))+a*b*t)/a

def sigma_tilde(t):
    return np.sqrt(sigma**2/(2*a**3)*(2*a*t + 4*np.exp(-a*t) - np.exp(-2*a*t) - 3))

def f_tau(t):
    return lambda_p * np.exp(-lambda_p * t)

def Integrand(t,z,csi,x,sign):
    """
    INPUT:
        t = row vector
        z = column vector
        csi = 
        x = 
    """
    z = z.reshape((len(z),1)) # z must be a column vector
    t = t[sigma_tilde(t) > 0]
    I = norm.pdf(np.log((z + sign * csi) / x),loc=mu_tilde(t),
     scale=sigma_tilde(t)) * f_tau(t)
    return I



    
def pfDESext2(z,x,u,J_jump):
    f = np.zeros(z.shape)
    csi = x*J_jump*u
    idx1 = z >= x+csi
    idx2 = z >= z-csi
    eta = 0.001
    t = np.arange(-1.6,2,eta)
    if z[idx1].size != 0:
        psi=lambda t: Integrand(np.exp(t - np.exp(-t)),z[idx1],csi,x,-1)*\
        np.exp(t - np.exp(-t))*(1+np.exp(-t))
        f[idx1] = p/(z[idx1]-csi)*integrate.trapz(dx=eta,y=psi(t))
          
    if z[idx2].size != 0:
        psi=lambda t: Integrand(np.exp(t - np.exp(-t)),z[idx2],csi,x,1)*\
        np.exp(t - np.exp(-t))*(1+np.exp(-t))
        f[idx2] = f[idx2]+(1-p)/(z[idx2]+csi)*integrate.trapz(dx=eta,y=psi(t))
    return f


u = 0.8
x = 1
eta = 0.0001
a = 0.5;b = 1.9
z = np.arange(a,b,eta)
t1 = time.time()
v = pfDESext2(z,x,u,J_jump)
t2 = time.time()
print(t2-t1)
plt.plot(z,v,'.b-')
I=integrate.simps(x = z,y = v)
print(I)