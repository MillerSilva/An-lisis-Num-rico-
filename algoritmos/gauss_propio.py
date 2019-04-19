import numpy as np


A=np.array([[2.0,-1.0,0.0],[1.0,6.0,-2.0],[4.0,-3.0,8.0]]).reshape(3,3)
b=np.reshape(np.array([2.0,-4.0,5.0]),(3,1))
n=np.shape(A)[0]
print("Matriz inicial: ")
print(A)
print("***********")
for k in range(0,n-1):
	for i in range(k+1,n):
		if A[i,k] != 0.0 :
			lam = A[i,k]/A[k,k]
			print("lam= "+str(lam))
			A[i,k:n] = A[i,k:n] - lam*A[k,k:n]
			b[i,:] = b[i,:] - lam*b[k,:]
	print(A)
	
print("Matriz final : ")			
print(A)	
print(b)


for k in range(n-1,-1,-1):
	b[k]=(b[k] - A[k,k+1:n].dot(b[k+1:n]) )/A[k,k]

print("Solucion final: ")
print(b)
