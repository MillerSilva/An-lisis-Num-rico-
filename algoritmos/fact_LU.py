import numpy as np

def suma(L,U,k,j):
    ac = 0
    for s in range(0,k):   # tener cuidado con este bucle, con el valor de k
        ac=ac+L[k,s]*U[s,j]
    return ac
    
    
def LU_doolitle(A):   
    """
        Metodo L1U sin hacer pivoteo
    """
    n=np.shape(A)[0]
    L=np.zeros((n,n))
    U=np.zeros((n,n)) 
    for k in range(0,n):
        L[k,k]=1
        for j in range(k,n):
            U[k,j]=A[k,j]-suma(L,U,k,j)
        for i in range(k+1,n):
            L[i,k]=(A[i,k]-suma(L,U,i,k))/U[k,k]
    print(L)
    print("******")
    print(U)                
    print("------")
    return L,U
    
def LU_crout(A):
    """
       Metodo LU1 sin hacer pivoteo
    """
    n=np.shape(A)[0]
    L=np.zeros((n,n))
    U=np.zeros((n,n)) 
    for k in range(0,n):
        U[k,k]=1
        for j in range(k,n):
            L[j,k]=A[j,k]-sum(L[j,p]*U[p,k] for p in range(k))
        for i in range(k+1,n):
            U[k,i]=(A[k,i]-sum(L[k,p]*U[p,i] for p in range(k)))/L[k,k]
    print(L)
    print("******")
    print(U)                
    print("------")    
    return L,U
    
### resolviendo LZ=b
def Lzb(L,b):
    n=np.shape(L)[0]
    z=np.zeros([n,1])
    z[0,0]=b[0,0]
    for i in range(1,n):
        z[i,0]= b[i,0] - sum(L[i,j]*z[j,0] for j in range(i))
        
    return z

### Ux=z

def Uxz(U,z):
    n=np.shape(U)[0]
    x=np.zeros([n,1])
    x[n-1,0]=z[n-1,0]/U[n-1,n-1]
    for i in range(n-2,-1,-1):
        x[i,0] = (z[i,0] - sum(U[i,j]*x[j,0] for j in range(i+1,n)))/U[i,i]
    return x

"""    
def LUsol(L,U,b):
    x=b
    n = len(x)
    # Solution of [L]{y} = {b}
    for k in range(n):
        x[k] = (b[k] - np.dot(L[k,0:k],b[0:k]))/L[k,k]
    # Solution of [L_transpose]{x} = {y}
    for k in range(n-1,-1,-1):
        x[k] = (b[k] - np.dot(U[k+1:n,k],b[k+1:n]))/U[k,k]
    return x    
    
    """

#A=np.reshape(np.array([10.0,10.0,20.0,20.0,25.0,40.0,30.0,50.0,61.0],dtype=float),(3,3))
A=np.reshape(np.array([4.0,-2.0,1.0,20.0,-7.0,12.0,-8.0,13.0,17.0],dtype=float),(3,3))
b=np.reshape(np.array([11.0,70.0,17.0]),(3,1))   
L1,U1=LU_doolitle(A)
L2,U2=LU_crout(A)    
print("Crout")
print(L2.dot(U2))
z = Lzb(L1,b)
print("z= ")
print(z)
x= Uxz(U1,z)
print("x= ")          
print(x)            
#print("Metodo alternativo: ")
#x_alt=LUsol(L1,U1,b)    
#print(x_alt)
