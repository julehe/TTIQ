"""
Returns a Chebyshev differentiation matrix.
The code was translated to Python based on the MATLAB script "cheb.m" given
in the book "Spectral Methods in MATLAB (2000)" by Trefethen. The MATLAB code
is publicly available on the book's homepage.
"""

import numpy as np

def cheb(N):
    if N==0: 
        D=0
        x=1
        return(D,x)
    else: 
        x = np.cos(np.pi*np.arange(N+1)/N); 
        c = np.array([2])
        for i in range(0,N-1):
            c = np.concatenate((c, np.array([1])), axis=0)
        c = np.concatenate((c, np.array([2])), axis=0)
        for i in range(0,len(c)):
            c[i] = c[i]*(-1)**i
        c.shape = (N+1,1)
        X = np.conjugate(np.array([x,]*(N+1))).T
        dX = X-np.conjugate(X).T
        D = c.dot(np.conjugate((1/c)).T) / (dX + (np.identity(N+1)))
        D = D - np.diag(sum(np.conjugate(D).T)) 
        return(D)
       

