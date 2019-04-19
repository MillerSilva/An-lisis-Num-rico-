#!/usr/bin/python
import numpy as np

def suma(L,d,i,j):
    ac = 0.0
    for v in range(0,j):
        ac=ac+L[i,v]*d[v,v]*L[j,v]
    return ac

def suma2(L,d,j):
    ac = 0.0
    for v in range(0,j):
        ac=ac+(d[v,v]*(L[j,v]*L[j,v]))
        print("***")
        print(ac)
        print("---")
    return ac


def LDL(A):
    n=np.shape(A)[0]
    L=np.zeros([n,n],dtype=float)
    d=np.zeros([n,n],dtype=float)
    for j in range(0,n):
        L[j,j] = 1
        d[j,j] = A[j,j] - suma2(L,d,j)
        print("Valor "+str(j+1)+" de la diagonal:")
        print(d[j,j])    
        for i in range(j+1,n):
            L[j,i] = 0
            L[i,j] = (A[i,j] - suma(L,d,i,j))/d[j,j]
    
    
    return L,d
    
    
A=np.reshape(np.array([4.0,3.0,2.0,1.0,3.0,3.0,2.0,1.0,2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0],dtype=float),(4,4))   
L,d=LDL(A)
print(L)
print("***************")
print(d)
print("Verifiquemos: ")
print((L.dot(d)).dot(L.T))
    
    
