import numpy as np
import scipy.linalg as sl
""""
def gauss_jordan(A): # retorna la inversa de la matriz A
    nrow = np.shape(A)[0] # numero de filas de A
    
    i= 1
    invA = np.identity(nrow)    #inicializa invA
    while i < nrow:
        a = np.array(A[:,i]/A[i,i]).reshape(nrow, 1)
        a[i] = 1-1/A[i,i]
        e = np.zeros((1, nrow))
        e[0,i] = 1
        
        L = np.identity(nrow) - np.dot(a, e)
        invA = np.dot(L, invA)
        if np.allclose(np.dot(invA, A),np.identity(nrow)):
            return invA
        i += 1
    print("La matriz :\n{}\n\nno es inversible".format(A))
"""

def max_norm(A):    # calcula la norma del maximo de la matriz A
    sum_row = sum(np.abs(A[0,:]))
    nrow = np.shape(A)[0]
    for i in range(1,nrow):
        s = sum(np.abs(A[i,:]))
        if s > sum_row:
            sum_row = s
    
    return sum_row


def cond(A):    # calcula la condicional de la matriz A
    invA = sl.inv(A) #DEPURAR GAUSS JORDAN , usando modulo scipy.linalg

    return max_norm(A)*max_norm(invA)

"""
tener cuidado con el tipo de datos del array(cantidades pequeñas lo aproxima a 0 
y asi para matrices no triangulares retorna True 
"""

# verifica si una matriz es triangular inferior teniendo en cuenta que los elementos bajo la diagonal difieren a lo sumo ERROR=1e-7 de cero
def istriangular_sup(A):   
    nrow = np.shape(A)[0]

    for k in range(0, nrow):
        if not np.allclose(A[k+1:, k] ,0.0) or  not np.isclose(A[nrow-1, nrow-1], 0.0):
            return False
    return True

def sup_triang(A):
    nrow=np.shape(A)[0]

    for i in range(nrow) :
        for j in range (nrow) :
            if i>j :
                if A[i,j]>1e-7 :
                    return False
    return True
# depurar
def eliminacion_gaussiana(A, b):
    nrow = np.shape(A)[0]
    
    i=0
    while not sup_triang(A) and i < nrow:
        #inicializa los vectores a y e
        a = np.zeros((nrow, 1))
        a[i+1:, i] = A[i+1:, i]/A[i,i]
        e = np.zeros((1, nrow))
        e[0, i] = 1

        L = np.identity(nrow) - np.dot(a, e)
        A = np.dot(L, A)
        b = np.dot(L, b)

        i+=1
    return A,b



def solve(A, b):
    #tA, tb: matrices obtenidas  despues de aplicar la transformacion gaussiana
    tA, tb = eliminacion_gaussiana(A, b)
    nrow = np.shape(A)[0]
    x = np.zeros(nrow)

    # resulve el sistema triangular superior (tA)x = tb, desde abajo hacia arriba
    for i in range(nrow-1, -1, -1): # range(nrow-1, -1, -1) = n-1, n-2, n-3, ..., 1, 0
        s = sum(tA[i,j]*x[j] for j in range(i+1, nrow))
        x[i] = (tb[i]- s)/tA[i,i] 
    
    return x

#depurar
def error_relativo(A, b):
    x0 = solve(A, b)    # es la aproximacion de la solucion sel sistema Ax=b, resuelto mediante eliminacion gaussiana
    r = np.dot(A, x0) - b   #vector residuo
    
    cota_sup = cond(A) * max_norm(r)/max_norm(b)
    cota_inf = max_norm(r)/(max_norm(b)*cond(A))

    print("El error relativo de la aproximacion({}) esta entre:\n{}\n{}".format(x0, cota_inf, cota_sup)) 

#depurar
def refinamiento_iterativos(A, b, MAX_ITERATIONS=10):
    x0 = solve(A, b)    # aproximacion de la solucion del sistema Ax=b
    invA =sl.inv(A)  #aproximacion de la inversa de A
    print("k | xk")
    print("---------------")
    print("{} | \t{}".format(0,x0))
    for k in range(MAX_ITERATIONS) :
        x=x0 + np.dot(invA,b-np.dot(A,x0))
        print("{} | \t{} ".format(k,x))
        x0=x