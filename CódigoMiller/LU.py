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


A=np.array([[8.04411657e-04, 1.66983660e-03, 4.17006102e-03, 3.92848931e-03,
        7.04645496e-03, 3.91197768e-03, 6.67480690e-03, 2.79508803e-03,
        2.22990408e-03, 6.97823920e-03, 4.33678963e-03, 7.43100852e-03,
        3.93985181e-03, 4.00524251e-03, 6.26144309e-03, 3.25852467e-04,
        3.15349399e-03, 5.04555795e-04, 1.35357205e-03, 6.99382935e-03],
       [2.04864259e-03, 7.27128167e-03, 6.41452566e-04, 6.19392264e-03,
        6.15321934e-03, 4.71903211e-04, 3.59093338e-03, 2.82185201e-03,
        8.56991596e-04, 1.34452803e-03, 3.35735583e-03, 3.68833039e-03,
        4.27109626e-03, 5.73995772e-03, 2.49194893e-03, 3.44445901e-03,
        4.09118764e-03, 2.03135769e-03, 6.79807534e-03, 5.07868401e-03],
       [5.87672919e-04, 4.38960775e-03, 4.60795943e-03, 2.23717545e-03,
        7.30269280e-03, 4.63825842e-03, 6.08457506e-03, 7.10584032e-03,
        2.55275656e-03, 6.42347089e-03, 7.16916712e-03, 6.92598616e-03,
        7.44320012e-03, 3.12082605e-03, 8.81388078e-04, 2.85653152e-03,
        3.74489197e-04, 6.26169138e-03, 3.83596291e-03, 1.42269702e-04],
       [7.98762564e-03, 7.05247819e-03, 7.47248135e-03, 1.67406519e-03,
        5.62945985e-03, 3.94921571e-03, 4.53273714e-03, 3.30875808e-03,
        5.15518898e-03, 7.42205011e-03, 5.06604511e-03, 7.28754521e-03,
        6.52293186e-03, 1.77064387e-03, 5.46471792e-03, 4.30901217e-03,
        4.93614544e-03, 6.26700266e-03, 2.93314995e-03, 5.99010293e-03],
       [5.91473213e-04, 7.61131303e-03, 1.82789290e-03, 2.89082149e-03,
        3.98217596e-03, 4.13656057e-03, 1.78640057e-03, 7.51838520e-03,
        3.41634853e-03, 4.31471198e-04, 2.12450084e-03, 3.03142521e-03,
        3.98349850e-03, 8.56551026e-04, 2.91179432e-03, 5.93272545e-03,
        7.26575960e-03, 3.96731122e-03, 3.09572290e-03, 4.74497206e-04],
       [4.07087708e-03, 7.27169982e-03, 2.92636972e-03, 1.86491580e-03,
        5.27902772e-03, 4.12774912e-03, 5.89166743e-03, 1.35691996e-03,
        6.80756348e-03, 7.49154906e-03, 9.75710399e-04, 4.61704451e-03,
        3.19674055e-03, 3.01969118e-03, 3.10737428e-03, 3.44276522e-03,
        1.83176082e-03, 2.99484542e-03, 4.75825617e-03, 2.34085984e-03],
       [7.88269778e-03, 8.59193026e-04, 5.28646381e-03, 1.21642647e-03,
        3.60231743e-03, 6.02518178e-03, 4.03540406e-03, 1.17240160e-05,
        5.63154712e-03, 1.54879435e-03, 7.32301084e-03, 3.16745141e-03,
        1.25157473e-04, 1.27444747e-03, 4.71873282e-03, 3.99888161e-04,
        2.41913480e-03, 2.70844522e-04, 1.95314374e-03, 5.70562060e-03],
       [7.76108413e-03, 2.17135844e-03, 1.31672500e-03, 3.81128931e-03,
        5.05658189e-03, 5.02004854e-03, 5.60031985e-03, 2.30625051e-04,
        3.84026569e-03, 3.71702981e-03, 2.57405504e-03, 6.05671972e-03,
        2.08940040e-03, 4.43829510e-03, 6.28396268e-03, 2.78031745e-03,
        4.89018859e-03, 2.61038865e-03, 6.99668765e-04, 7.27623472e-03],
       [7.71470791e-03, 6.19447797e-03, 3.24510055e-03, 4.85557406e-03,
        1.51927909e-03, 5.24749016e-03, 6.46264156e-03, 6.87182296e-03,
        3.27335153e-03, 1.42419893e-03, 3.00842631e-04, 5.11023644e-03,
        5.93095246e-04, 2.54443849e-03, 5.57369590e-03, 4.06959977e-03,
        7.03416004e-03, 4.62556134e-03, 6.49885563e-03, 7.09935353e-03],
       [2.51197502e-03, 7.62180113e-03, 7.25070873e-03, 2.33338078e-03,
        4.53225124e-04, 5.96562693e-04, 5.87562396e-03, 4.67156163e-03,
        3.77959154e-03, 3.70751209e-03, 4.23868969e-03, 4.92367750e-03,
        1.88722327e-03, 6.52834043e-03, 1.59152760e-03, 1.21099951e-03,
        3.54580409e-03, 2.66138580e-04, 6.25229062e-03, 5.22449257e-03],
       [5.05664027e-03, 1.77966725e-03, 6.97960304e-03, 6.68980538e-03,
        6.14040543e-03, 5.04558765e-03, 5.12662385e-03, 6.34416124e-03,
        3.23284987e-03, 5.69905847e-03, 1.53460688e-03, 3.08815830e-03,
        1.74363311e-03, 5.60092293e-03, 3.36801985e-03, 1.83599746e-03,
        3.67778794e-03, 6.96542953e-03, 3.01452537e-03, 9.93852002e-04],
       [6.99481564e-03, 6.31870740e-03, 4.29175044e-03, 7.51304941e-03,
        1.04953882e-03, 1.74854773e-03, 2.60100715e-03, 5.34540370e-03,
        4.57215146e-04, 7.00580753e-03, 7.27801217e-03, 3.97081116e-03,
        3.30436751e-03, 4.71941158e-03, 6.83768509e-03, 7.03341702e-03,
        6.65235069e-03, 3.10109100e-03, 4.14226337e-03, 2.72470794e-03],
       [4.72661427e-05, 6.19387130e-03, 6.58376701e-03, 1.65593049e-03,
        5.25684151e-03, 7.49742282e-03, 3.42107454e-03, 1.86505707e-05,
        4.94834131e-04, 3.54516496e-03, 9.71965406e-04, 5.41346795e-04,
        1.84427036e-03, 1.43861600e-03, 3.13265641e-03, 1.49243318e-03,
        9.75205910e-04, 2.89985852e-03, 3.57751036e-03, 3.87848181e-03],
       [5.49607939e-03, 5.67428139e-04, 3.04196553e-03, 2.15170719e-03,
        3.33220478e-03, 4.31235201e-03, 5.33438865e-03, 6.04143372e-03,
        1.14554342e-03, 3.94838595e-03, 8.00733329e-04, 2.93666824e-03,
        5.91244347e-04, 6.69993634e-03, 1.83499351e-03, 2.69668640e-03,
        1.71102665e-03, 5.40636188e-03, 1.98107807e-03, 5.04972417e-03],
       [7.23276661e-03, 3.75844577e-03, 5.22378441e-03, 8.11000633e-04,
        4.42875119e-03, 7.69082627e-03, 2.40981558e-03, 2.44630678e-03,
        8.61007357e-04, 5.23112688e-03, 5.92412372e-03, 1.39220873e-03,
        2.76379162e-03, 1.98805217e-03, 4.33564783e-03, 5.78389770e-03,
        1.61014446e-04, 1.48978125e-03, 4.07568457e-03, 4.19326104e-04],
       [3.28260105e-03, 5.99409966e-03, 5.18988066e-03, 3.35423772e-03,
        7.59196028e-04, 1.11066121e-04, 6.80048074e-03, 3.61619199e-04,
        7.15321783e-03, 3.96939291e-03, 1.87027178e-03, 6.80150467e-03,
        6.76324459e-03, 6.93119573e-03, 3.47782983e-04, 4.16203769e-03,
        6.44205081e-03, 7.19447414e-03, 9.28206898e-04, 1.51511425e-03],
       [3.07262906e-03, 3.57279499e-03, 4.39519796e-04, 2.98994261e-03,
        6.90437203e-03, 2.18452527e-03, 7.05915239e-04, 6.91730450e-03,
        6.13915514e-03, 7.51487981e-03, 1.46222184e-03, 7.71954340e-04,
        4.61142735e-03, 2.70000132e-03, 3.70079237e-03, 1.85291549e-03,
        6.28674341e-03, 1.06459613e-03, 4.97308272e-04, 6.25541164e-03],
       [5.37956111e-03, 3.30958767e-03, 2.23824409e-04, 2.02443450e-03,
        3.51203282e-03, 2.80452331e-03, 6.79466147e-03, 1.93514120e-03,
        4.06745127e-03, 8.20126182e-04, 1.54810444e-03, 6.75925668e-03,
        3.42020451e-03, 5.04713078e-03, 3.31264385e-03, 4.21679619e-03,
        5.04612586e-03, 7.34077667e-03, 5.43190169e-03, 1.55053720e-03],
       [5.46731936e-03, 1.94436872e-03, 4.02267391e-03, 3.24190006e-03,
        2.24988294e-03, 1.58750990e-03, 3.65624478e-03, 7.40794844e-03,
        3.96425627e-03, 5.56059572e-03, 4.57399018e-03, 2.98205792e-03,
        1.46649334e-03, 1.20829948e-03, 7.13753184e-03, 2.89452782e-04,
        3.52344962e-03, 4.79919747e-03, 1.93149038e-03, 1.05995127e-03],
       [1.71119989e-04, 3.53328733e-03, 3.55660851e-03, 4.75816056e-03,
        4.87744370e-03, 4.14596580e-03, 3.43022364e-03, 5.34508446e-03,
        4.20477771e-03, 5.08690714e-03, 5.95031043e-03, 6.37161536e-03,
        6.57471871e-03, 3.57066130e-03, 4.75554427e-03, 1.96297083e-03,
        3.61848417e-03, 3.04132691e-03, 4.27280177e-03, 9.88100101e-04]])
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