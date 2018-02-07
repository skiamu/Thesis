# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:38:54 2018

@author: andrea
"""

from pyOpt import Optimization
from pyOpt import SLSQP
import numpy as np


def objfun(x,**kwargs):
    W = kwargs['W']
    a = kwargs['a']
    f = a * W.dot(x).dot(x)
    g = [x.sum()-1]
    fail = 0
    return f,g,fail
    
    
opt_prob = Optimization('test problem',objfun)
opt_prob.addVarGroup('x',3,'c',lower=np.zeros(3),
                     upper=np.ones(3),value=[1,0,0])
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
print opt_prob


W = np.eye(3)


# Instantiate Optimizer (SLSQP) & Solve Problem
slsqp = SLSQP()
slsqp.setOption('IPRINT',-1)
slsqp(opt_prob,sens_type='FD', W=W,a=10)
print opt_prob.solution(0)








    
    
    