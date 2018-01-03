#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 11:42:01 2017

@author: andrea
"""
import numpy as np

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
    return (D,DiscreteS)      
            
    
def objfun(x,arg):
    """
    INPUT:
        x = vector of parameter
        arg = (D,dt)
    """
    D = arg[0]
    dt = arg[1]
    sign = arg[2] # 1 for minimization, -1 for maximization
    Nzero = len(D[D==0])
    Nuno = len(D[D==1])
    NmenoUno = len(D[D==-1])
    f = -Nzero*x[1]*dt + Nuno*np.log((1-np.exp(-x[1]*dt))*x[0])+ \
        NmenoUno*np.log((1-np.exp(-x[1]*dt))*(1-x[0]))
    return sign * f
    
    
    
    
    
    