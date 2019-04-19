import numpy as np

"""
El programa debe de recibir la matriz A triangular superior y el vector b con
todas las operaciones realizadas sobre el
"""
def solve_triinf(A,b)
	x=b
	for k in range(n-1,-1,-1):
		x[k]=(x[k] - A[k,k+1:n].dot(x[k+1:n]) )/A[k,k]
	return x



	
