import numpy as np

def swap(A,r1,r2):
	A[[r1,r2]]=	A[[r2,r1]]

def suma(L,U,k,j):
    ac = 0
    for s in range(0,k):        # tener cuidado con este bucle, con el valor de k
        ac=ac+L[k,s]*U[s,j]
    return ac
    
def posee_LU(A):
    n=np.shape(A)[0]    
    for k in range(0,n):
        if np.linalg.det(A[k:,k:]) == 0.0:
            return False
    return True
    
def LU_doolitle_piv(A):  
    if posee_LU(A):
        """
            Metodo L1U con pivoteo
        """
        n=np.shape(A)[0]
        L=np.zeros((n,n))
        U=np.zeros((n,n)) 
        for k in range(0,n):
            L[k,k]=1
            for j in range(k,n):
                U[k,j]=A[k,j]-suma(L,U,k,j)
            
            #p=np.argmax(abs(A[k:n,k]))+k
            if p != k:
                print("Pivoteando")
                swap(A,p,k)
                
            for i in range(k+1,n):
                L[i,k]=(A[i,k]-suma(L,U,i,k))/U[k,k]
        print(L)
        print("******")
        print(U)                
        print("------")
        return L,U
    else:
        print("No posee LU")
    
def LU_crout_piv(A):
    if posee_LU(A):    
        """
           Metodo LU1 con pivoteo
        """
        n=np.shape(A)[0]
        L=np.zeros((n,n))
        U=np.zeros((n,n)) 
        for k in range(0,n):
            U[k,k]=1
            for j in range(k,n):
                L[j,k]=A[j,k]-suma(U,L,k,j)
                
            p=np.argmax(abs(A[k:n,k]))+k
            if p != k:
                #print("Pivoteando")
                swap(A,p,k)
            
            for i in range(k+1,n):
                U[k,i]=(A[k,i]-suma(L,U,k,i))/L[k,k]
        print(L)
        print("******")
        print(U)                
        print("------")    
        return L,U

    else:
        print("No posee LU")
    
def LUsolve(a,b):
    n = len(a)
    for k in range(1,n):
        b[k] = b[k] - np.dot(a[k,0:k],b[0:k])

    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    return b


    
A=np.reshape(np.array([4.0,-8.0,1.0,-6.0,4.0,-2.0,1.0,-2.0,3.0]),(3,3))
L1,U1=LU_doolitle_piv(A)
L2,U2=LU_crout_piv(A)    
                
