"""
Implementacion del metodo LU
"""

import numpy as np
import scipy.linalg as sl

# retorna el primer indice del elemento mayor(si hubiera dos mayores)
def max_index(x):   
    element_max = max(x)
    index = [k for k in range(len(x)) if x[k] == element_max]
    return index[0]

# cambia la fila k por la fila j en A
def change_row(nrow,k, j):
    A=np.identity(nrow)
    A[k,] = A[j,]
    A[j,] = np.identity(nrow)[k,]
    return A

#Para resolver sistema triangular superior

def sumS(k,nrow,A,x):
    n=nrow-1
    sums=0
    if k==0:
        return sums
    else:
        for r in range(n-k+1,nrow):
            sums=sums+A[n-k][r]*x[r][0]
        return sums

def SupSolve(A,b):
    nrow=np.shape(A)[0]
    n=nrow-1
    x=np.zeros(nrow).reshape(nrow, 1)
    for k in range(nrow):
        x[n-k][0]=(b[n-k][0]-sumS(k,nrow,A,x))/A[n-k][n-k]
    return x    
    SupSolve(A,b)
    
    
#Para resolver sistema triangular inferior    

def sumI(k,A,x):
    sumi=0
    if k==0:
        return sumi
    else:
        for r in range(k):
            sumi=sumi+A[k][r]*x[r][0]
        return sumi 

def InfSolve(A,b):
    nrow=np.shape(A)[0]
    x=np.zeros(nrow).reshape(nrow, 1)
    for k in range(nrow):
        x[k][0]=(b[k][0]-sumI(k,A,x))/A[k][k]
    return x

#retorna la factorizacion LU de A 
def lu(A):
        nrow = np.shape(A)[0]
        inverL = np.identity(nrow)
        U = A
        P = np.identity(nrow)
        for k in range(nrow-1):
            index = max_index(U[k:,k])+k
            #print(index,k)
            # crea P (matrix de permutacion), iteracion k
            Pk = change_row(nrow, k, index)
            #print(Pk)
            P = np.dot(Pk, P)
            #
            U=np.dot(Pk,U)
            #print(U)
            # crea L (elimina elementos en la fila k), iteracion k
            a = np.zeros(nrow).reshape(nrow, 1)
            a[k+1:,0] = U[k+1:,k]/U[k,k]
            e = np.zeros(nrow).reshape(1, nrow)
            e[0,k] = 1
            Lk = np.identity(nrow) - np.dot(a, e)
            #print(Lk)
            U = np.dot(Lk, U)
            #print(U)
            LPk = np.dot(Lk,Pk)
            inverL = np.dot(LPk, inverL)
        L = np.dot(P, sl.inv(inverL)) #calcula la inversa de una matriz triangular (inverL)

        return L, U, P
    
# resuelve el sistema LUx = b
def solvelu(A, b):
    if sl.det(A)!=0:
        L, U, P = lu(A)
        z = InfSolve(L,np.dot(P,b)) # resuelve el sistema Lz = b, donde z = Ux y L es triangular inferior

        x = SupSolve(U, z)  # resuelve el sistema Ux = z,  donde U es triangular superior

        return x
    else:
        print("La matriz no tiene descomposicion LU.")
        return np.zeros_like(b)    
    
B=np.array([1,1,1,1,2,3,2,7,3,1,7,1,10,11,12,1]).reshape(4,4)
c=np.array([20,80,60,154]).reshape(4,1)
solvelu(B,c)
C=np.array([1,1,1,2,1,4,2,5,6]).reshape(3,3)
d1=np.array([6,16,30]).reshape(3,1)
#solvelu(C,d1)


A=np.array([[8.044e-04, 1.669e-03, 4.170e-03, 3.928e-03, 7.046e-03, 3.911e-03, 6.674e-03, 2.795e-03, 2.229e-03,6.978e-03, 
	4.336e-03, 7.431e-03, 3.939e-03, 4.005e-03, 6.261e-03, 3.258e-04, 3.153e-03, 5.045e-04, 1.353e-03, 6.993e-03],
       [2.048e-03, 7.271e-03, 6.414e-04, 6.193e-03, 6.153e-03, 4.719e-04, 3.590e-03, 2.821e-03, 8.569e-04, 1.344e-03, 
	3.357e-03, 3.688e-03, 4.271e-03, 5.739e-03, 2.491e-03, 3.444e-03, 4.091e-03, 2.031e-03, 6.798e-03, 5.078e-03],
       [5.876e-04, 4.389e-03, 4.607e-03, 2.237e-03, 7.302e-03, 4.638e-03, 6.084e-03, 7.105e-03, 2.552e-03, 6.423e-03,
 	7.169e-03, 6.925e-03, 7.443e-03, 3.120e-03, 8.813e-04, 2.856e-03, 3.744e-04, 6.261e-03, 3.835e-03, 1.422e-04],
       [7.987e-03, 7.052e-03, 7.472e-03, 1.674e-03, 5.629e-03, 3.949e-03, 4.532e-03, 3.308e-03, 5.155e-03, 7.422e-03, 
	5.066e-03, 7.287e-03, 6.522e-03, 1.770e-03, 5.464e-03, 4.309e-03, 4.936e-03, 6.267e-03, 2.933e-03, 5.990e-03],
       [5.914e-04, 7.611e-03, 1.827e-03, 2.890e-03, 3.982e-03, 4.136e-03, 1.786e-03, 7.518e-03, 3.416e-03, 4.314e-04, 
	2.124e-03, 3.031e-03, 3.983e-03, 8.565e-04, 2.911e-03, 5.932e-03, 7.265e-03, 3.967e-03, 3.095e-03, 4.744e-04],
       [4.070e-03, 7.271e-03, 2.926e-03, 1.864e-03, 5.279e-03, 4.127e-03, 5.891e-03, 1.356e-03, 6.807e-03, 7.491e-03, 
	9.757e-04, 4.617e-03, 3.196e-03, 3.019e-03, 3.107e-03, 3.442e-03, 1.831e-03, 2.994e-03, 4.758e-03, 2.340e-03],
       [7.882e-03, 8.591e-04, 5.286e-03, 1.216e-03, 3.602e-03, 6.025e-03, 4.035e-03, 1.172e-05, 5.631e-03, 1.548e-03, 
	7.323e-03, 3.167e-03, 1.251e-04, 1.274e-03, 4.718e-03, 3.998e-04, 2.419e-03, 2.708e-04, 1.953e-03, 5.705e-03],
       [7.761e-03, 2.171e-03, 1.316e-03, 3.811e-03, 5.056e-03, 5.020e-03, 5.600e-03, 2.306e-04, 3.840e-03, 3.717e-03, 
	2.574e-03, 6.056e-03, 2.089e-03, 4.438e-03, 6.283e-03, 2.780e-03, 4.890e-03, 2.610e-03, 6.996e-04, 7.276e-03],
       [7.714e-03, 6.194e-03, 3.245e-03, 4.855e-03, 1.519e-03, 5.247e-03, 6.462e-03, 6.871e-03, 3.273e-03, 1.424e-03, 
	3.008e-04, 5.110e-03, 5.930e-04, 2.544e-03, 5.573e-03, 4.069e-03, 7.034e-03, 4.625e-03, 6.498e-03, 7.099e-03],
       [2.511e-03, 7.621e-03, 7.250e-03, 2.333e-03, 4.532e-04, 5.965e-04, 5.875e-03, 4.671e-03, 3.779e-03, 3.707e-03, 
	4.238e-03, 4.923e-03, 1.887e-03, 6.528e-03, 1.591e-03, 1.210e-03, 3.545e-03, 2.661e-04, 6.252e-03, 5.224e-03],
       [5.056e-03, 1.779e-03, 6.979e-03, 6.689e-03, 6.140e-03, 5.045e-03, 5.126e-03, 6.344e-03, 3.232e-03, 5.699e-03, 
	1.534e-03, 3.088e-03, 1.743e-03, 5.600e-03, 3.368e-03, 1.835e-03, 3.677e-03, 6.965e-03, 3.014e-03, 9.938e-04],
       [6.994e-03, 6.318e-03, 4.291e-03, 7.513e-03, 1.049e-03, 1.748e-03, 2.601e-03, 5.345e-03, 4.572e-04, 7.005e-03, 
	7.278e-03, 3.970e-03, 3.304e-03, 4.719e-03, 6.837e-03, 7.033e-03, 6.652e-03, 3.101e-03, 4.142e-03, 2.724e-03],
       [4.726e-05, 6.193e-03, 6.583e-03, 1.655e-03, 5.256e-03, 7.497e-03, 3.421e-03, 1.865e-05, 4.948e-04, 3.545e-03, 
	9.719e-04, 5.413e-04, 1.844e-03, 1.438e-03, 3.132e-03, 1.492e-03, 9.752e-04, 2.899e-03, 3.577e-03, 3.878e-03],
       [5.496e-03, 5.674e-04, 3.041e-03, 2.151e-03, 3.332e-03, 4.312e-03, 5.334e-03, 6.041e-03, 1.145e-03, 3.948e-03, 
	8.007e-04, 2.936e-03, 5.912e-04, 6.699e-03, 1.834e-03, 2.696e-03, 1.711e-03, 5.406e-03, 1.981e-03, 5.049e-03],
       [7.232e-03, 3.758e-03, 5.223e-03, 8.110e-04, 4.428e-03, 7.690e-03, 2.409e-03, 2.446e-03, 8.610e-04, 5.231e-03, 
	5.924e-03, 1.392e-03, 2.763e-03, 1.988e-03, 4.335e-03, 5.783e-03, 1.610e-04, 1.489e-03, 4.075e-03, 4.193e-04],
       [3.282e-03, 5.994e-03, 5.189e-03, 3.354e-03, 7.591e-04, 1.110e-04, 6.800e-03, 3.616e-04, 7.153e-03, 3.969e-03, 
	1.870e-03, 6.801e-03, 6.763e-03, 6.931e-03, 3.477e-04, 4.162e-03, 6.442e-03, 7.194e-03, 9.282e-04, 1.515e-03],
       [3.072e-03, 3.572e-03, 4.395e-04, 2.989e-03, 6.904e-03, 2.184e-03, 7.059e-04, 6.917e-03, 6.139e-03, 7.514e-03, 
	1.462e-03, 7.719e-04, 4.611e-03, 2.700e-03, 3.700e-03, 1.852e-03, 6.286e-03, 1.064e-03, 4.973e-04, 6.255e-03],
       [5.379e-03, 3.309e-03, 2.238e-04, 2.024e-03, 3.512e-03, 2.804e-03, 6.794e-03, 1.935e-03, 4.067e-03, 8.201e-04, 
	1.548e-03, 6.759e-03, 3.420e-03, 5.047e-03, 3.312e-03, 4.216e-03, 5.046e-03, 7.340e-03, 5.431e-03, 1.550e-03],
       [5.467e-03, 1.944e-03, 4.022e-03, 3.241e-03, 2.249e-03, 1.587e-03, 3.656e-03, 7.407e-03, 3.964e-03, 5.560e-03, 
	4.573e-03, 2.982e-03, 1.466e-03, 1.208e-03, 7.137e-03, 2.894e-04, 3.523e-03, 4.799e-03, 1.931e-03, 1.059e-03],
       [1.711e-04, 3.533e-03, 3.556e-03, 4.758e-03, 4.877e-03, 4.145e-03, 3.430e-03, 5.345e-03, 4.204e-03, 5.086e-03, 
	5.950e-03, 6.371e-03, 6.574e-03, 3.570e-03, 4.755e-03, 1.962e-03, 3.618e-03, 3.041e-03, 4.272e-03, 9.881e-04]])
d=np.array([[10000.],
       [99000.],
       [36000.],
       [97000.],
       [45000.],
       [61000.],
       [82000.],
       [74000.],
       [34000.],
       [51000.],
       [79000.],
       [62000.],
       [42000.],
       [28000.],
       [71000.],
       [52000.],
       [58000.],
       [33000.],
       [31000.],
       [31000.]])       
I=np.identity(20)
x=solvelu(I-A,d)
dr=np.dot(I-A,x)
R=dr-d
dr
R
x
       
