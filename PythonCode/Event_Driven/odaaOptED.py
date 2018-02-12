# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:48:00 2018

@author: andrea
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 19:09:02 2018

@author: andrea
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 16:19:02 2018

@author: andrea
"""
import numpy as np
from scipy.interpolate import interp1d
from pyOpt import Optimization
from pyOpt import SLSQP # (*)
from pyOpt import PSQP
from pyOpt import CONMIN
from pyOpt import COBYLA
from pyOpt import SOLVOPT
from pyOpt import KSOPT
from pyOpt import NSGA2
from pyOpt import SDPEN
from pyOpt import KSOPT # (*)
from pyOpt import ALGENCAN
from pyOpt import SDPEN
import matplotlib.pyplot as plt

# ========================================================================
#   FUNCTIONs DEFINITION
# ========================================================================
   
def ODAAalgorithmED(N,X,model):
    ''' 
    INPUT:
        N = number of time steps
        X = discretized target sets [tuple dimension N+1]
        model = object containing model information and parameters
    OUTPUT:
        J = optial value function at each tie instant [list]
        U = optimal asset allocation [list of np array]
    '''
    U = [0]*N 
    J = [0]*(N+1)
    J[-1] = np.ones(len(X[-1]))
    eta = 1e-5/5 # discretization step
    for k in range(N-1,N-4,-1):
        print k # print current iteration
        u0 = -1
        dimXk = len(X[k])
        Uk = np.zeros(dimXk)
        Jk = np.zeros(dimXk)
        int_domain = np.arange(X[k+1][0],X[k+1][-1],eta)
        Jinterp = interp1d(X[k+1],J[k+1])(int_domain)
        for j in range(dimXk-1,-1,-1):
            print j
            opt_prob = solveOpt(int_domain,Jinterp,X[k][j],model,u0,-1)
            Uk[j] = opt_prob.solution(0).getVar(0).value
            Jk[j] = opt_prob.solution(0).getObj(0).value
            u0 = Uk[j]
        U[k] = Uk
        J[k] = -Jk
        plt.figure()
        plt.stackplot(X[k],Uk)
        plt.title('k = '+ str(k))
    return J,U
    
    
def solveOpt(int_domain,J,a,model,u0,sign):  
    '''
    INPUT:
        int_domain = 
        J = 
        a = 
        model = 
        u0 =
        sign = 
    OUTPUT:
        opt_prob = 
    '''         
    def objfun(u,**kwargs):
        '''objfun defines optimization problem using the pyOpt sintax'''
        # 1) extract paraeters
        int_domain = kwargs['int_domain'] 
        J = kwargs['J'] 
        x = kwargs['a'] 
        model = kwargs['model'] 
        sign = kwargs['sign'] 
        # 2) define objective function and constraints
        funz = np.trapz(x = int_domain,y = J * model.pf(int_domain,u,x))
        g = []
        fail = 0
        return sign * funz,g,fail
    opt_prob = Optimization('ODAA problem',objfun)
    opt_prob.addObj('funz')
    solver = SLSQP() # choose the solver
    solver.setOption('IPRINT',0)
    opt_prob.addVar('u','c',lower=-1,upper=1,value=u0)
    #print opt_prob # print optimization problem                  
    solver(opt_prob,sens_type='FD',int_domain=int_domain,J=J,
          a=a,model=model,sign = sign)
    #print opt_prob.solution(0) # print solution
    return opt_prob
      
    
    