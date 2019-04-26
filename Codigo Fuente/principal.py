# -*- coding: utf-8 -*-
from toolNick import *      # Se encuentran herramientos para los metodos
from metodosNumericos import *   # Se encuentran todos los metodos
import numpy as np
from time import time       #Para el tiempo de ejecucion

def get_Data(A,b):
    """Esta funcion retorna un dicionario de datos como para cada diccionario de metodos con y sin pivoteo
    name_data = ['Algoritmo', 'Solucion','Error','Time']
   
    name_method = ['Eliminacion Gaussian','Gauss Jordan','Eliminacion LU','Grout',
                    'Doolittle','Factorizacion LDLt', 'Factorizacion Cholesky',
                    'Eliminacion Gaussian Pivote','Gauss Jordan Pivoteo','Eliminacion LU Pivoteo'
                    ]
    diccionary['name_method']['name_data']"""
    
    algoritmo = [elimGauss, gauss_Jordan, elim_LU, fac_Grout_LU1, fac_Grout_L1U, LDLt, cholesky,
                elimGauss, gauss_Jordan, elim_LU]
    name_method = ['Eliminacion Gaussian','Gauss Jordan','Eliminacion LU','Grout','Doolittle', 
                   'Factorizacion LDLt', 'Factorizacion Cholesky',
                   'Eliminacion Gaussian Pivote','Gauss Jordan Pivoteo','Eliminacion LU Pivoteo']
    fil,col = A.shape
    #Creamos los diccionarios a partir de lista de nombres
    method = {}
    cant_Method = len(name_method)
    i = 0
    for nameM in name_method:
        data = {}
        start_time = time()
        # Actualizamos el diccionario dato
        data.update({'Algoritmo' : algoritmo[i]})
        if i < cant_Method - 3:    # Sin pivote
            data.update({'Solucion' : data['Algoritmo'](A, b)})
        else:                     # Con pivote
            data.update({'Solucion' : data['Algoritmo'](A, b, pp=True)})	
        data.update({'Error': error_relativo(A,b.reshape(fil,1),data['Solucion'].reshape(fil,1))})
        time_end = time() - start_time
        data.update({'Time': time_end})
        #data.up
        # Actualizamos el diccionario method
        method.update({nameM : data})
        i += 1


    return method

if __name__ == '__main__' :
    A = np.array([[5, 1, -2, 0],
              [1, 8, 2, 0],
              [-2, 2, 7, 1],
              [0, 0, 1, 3]])

    b = np.array([10, 23, 5, 17])
    datos = get_Data(A,b)
    print(datos)
    help(get_Data)