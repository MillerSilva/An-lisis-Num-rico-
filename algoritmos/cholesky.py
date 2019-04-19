#!/usr/bin/python
import numpy as np
import math

def suma1(L,k):
    ac=0.0
    for s in range(0,k-1):
        ac=ac+(L[k,s]*L[k,s])
    return ac
        
def suma2(L,i,k):
    ac=0.0
    for s in range(0,k-1):
        ac=ac+(L[i,s]*L[k,s])    
    return ac

def cholesky_propio(A):
    n=np.shape(A)[0]
    L=np.zeros((n,n))
    for k in range(0,n):
       
       
        L[k,k]=(A[k,k]-sum(L[k][j] * L[k][j] for j in range(k)))**0.5
        
        for i in range(k+1,n):
           
            L[i,k]=( A[i,k] -  sum(L[i,s]*L[k,s] for s in range(0,k)) ) /L[k,k] 

            
    return L

def cholesky(A):
    """Performs a Cholesky decomposition of A, which must 
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix, L."""
    n = len(A)

    # Create zero matrix for L
    L=np.zeros((n,n))

    # Perform the Cholesky decomposition
    for i in range(n):
        for k in range(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in xrange(k))
            
            if (i == k): # Diagonal elements

                L[i][k] = math.sqrt(A[i][i] - tmp_sum)
            else:

                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L


def choleskiSol(L,b):
    x=b
    n = len(x)
    # Solution of [L]{y} = {b}
    for k in range(n):
        x[k] = (b[k] - np.dot(L[k,0:k],b[0:k]))/L[k,k]
    # Solution of [L_transpose]{x} = {y}
    for k in range(n-1,-1,-1):
        x[k] = (b[k] - np.dot(L[k+1:n,k],b[k+1:n]))/L[k,k]
    return x





A=np.reshape(np.array([1.44,-0.36,5.52,0.0,-0.36,10.33,-7.78,0.0,5.52,-7.78,28.40,9.0,0.0,0.0,9.0,61.0],dtype=float),(4,4))
b=np.reshape(np.array([0.04,-2.15,0.0,0.88]),(4,1))  
print(cholesky_propio(A))
L=cholesky_propio(A)
print(choleskiSol(L,b))
 
