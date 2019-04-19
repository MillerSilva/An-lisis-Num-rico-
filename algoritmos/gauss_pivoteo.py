# -*- coding: utf-8 -*-
import numpy as np

# FUNCION PARA INTERCAMBIAR FILAS
def swap(A,r1,r2):
	A[[r1,r2]]=	A[[r2,r1]]
	
A=np.reshape(np.array([4.0,-2.0,1.0,-6.0,4.0,-2.0,1.0,-2.0,3.0]),(3,3))
b=np.reshape(np.array([11.0,-16.0,17.0]),(3,1))
# Valor de tolerancia, esto se usara para realizar el pivoteo
tol = 1.0e-7 


n=np.shape(A)[0]

# arreglo que contendra los valores mas grandes de cada fila
s=np.zeros(n)

# aqui se rellena el arreglo s
for i in range(0,n):
    s[i] = max(np.abs(A[i,:]))
    
print("s= ")
print(s)


for k in range(0,n-1):
    
    p = np.argmax(np.abs(A[k:n,k])) + k
    if abs(A[p,k]) < tol:
        print("Error: matriz singular")
        break
    if p!=k:
        swap(A,k,p)
        swap(b,k,p)
    
    #Eliminacion    
    for i in range(k+1,n):
        if A[i,k] != 0.0 :
		    lam = A[i,k]/A[k,k]
		#print("lam= "+str(lam))
        A[i,k:n] = A[i,k:n] - lam*A[k,k:n]
        b[i,:] = b[i,:] - lam*b[k,:]
	
if abs(A[n-1,n-1]) < tol: 
    print("Error: matriz singular")
    break

print("******------")
print(A)
print("------*******")
print(b)			


