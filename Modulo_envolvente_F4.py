# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:53:21 2018

@author: Team REOS

Módulo para sacar la envolvente de vuelo del avión F4. Como input tiene el 
PESO del avión cargado con el REOS. Se realiza una interpolación con las
envolventes conocidas de unos pesos determinados y se realiza una aproximación
polinómica de alto orden para recrear las ecuaciones de los límites de la 
envolvente para el PESO de entrada.
"""

#from inputs_iniciales import MASS
from numpy import transpose, zeros, matmul, array, size, dot
from numpy.linalg import cholesky, inv

# Definción de variables. Todos los pesos en libras
#PESO = MASS # Para cuando quitemos la línea de abajo e importemos directamente del otro módulo
PESO = 47300.0  # Input de entrada. Peso del avión cargado con el cohete REOS
P1 = 42777.0 
P2 = 43035.0
P3 = 45472.0
P4 = 46279.0
P5 = 46537.0
P6 = 48974.0

alt_enve = []  # Lista de valores de altitud de la curva límite de la envolvente
mach_enve = []  # Lista de valores de Mach de la curva límite de la envolvente

F = open("Envolventes.txt", 'r')  # Apertura del fichero de lectura.

for i in F:  # Bucle para recorrer el archivo de texto
    i = i.strip()  # Sentencias para eliminar los espacios sobrantes de F
    i = i.split() 
    for j, a in enumerate(i):  # Bucle para recorrer las filas del archivo
            i[j] = float(i[j])  # Convierte los valores de strings a floats
    # Condicionales para separar los casos, las curvas entre las que hay que
    # interpolar. Interpolaciones lineales en altitud para mismo Mach
    # Para el último punto se realiza una interpolación lineal en el número
    # de Mach para igual valor final de altitud (36000ft)
    # Para PESO < P1 y para PESO > P6 lo que se realiza es una extrapolación 
    
    if PESO < P1:
        if i[0] < 0.9:
            cte = -0.000478303
        elif 0.9 <= i[0] >= 1:
            cte = -0.000322528
        elif 1 < i[0]:
            cte = -0.001747212
        
        enve = i[1] + cte * (PESO - P1)
        alt_enve.append(enve)
        mach_enve.append(i[0])

    elif P1 < PESO < P2:
        enve = i[1] + (i[2]-i[1]) * (PESO - P1) / (P2 - P1)
        alt_enve.append(enve)
        mach_enve.append(i[0])

    elif P2 < PESO < P3:
        if i[0] <= 1.88:
            enve = i[2] + (i[3] - i[2]) * (PESO - P2) / (P3 - P2)
            alt_enve.append(enve)
            mach_enve.append(i[0])                            
        else:
            Mfin = i[8] + (i[9] - i[8]) * (PESO - P2) / (P3 - P2)            
    
    elif P3 < PESO < P4:
        if i[0] <= 1.78:
            enve = i[3] + (i[4] - i[3]) * (PESO - P3) / (P4 - P3)
            alt_enve.append(enve)
            mach_enve.append(i[0])          
        else:
            Mfin = i[9] + (i[10] - i[9]) * (PESO - P3) / (P4 - P3)

    elif P4 < PESO < P5:
        if i[0] <= 1.76:
            enve = i[4] + (i[5] - i[4]) * (PESO - P4) / (P5 - P4)
            alt_enve.append(enve)
            mach_enve.append(i[0])
        else:
            Mfin = i[10] + (i[11] - i[10]) * (PESO - P4) / (P5 - P4)

    elif P5 < PESO < P6:
        if i[0] <= 1.5:
            enve = i[5] + (i[6] - i[5]) * (PESO - P5) / (P6 - P5)
            alt_enve.append(enve)
            mach_enve.append(i[0])
        else:
            Mfin = i[11] + (i[12] - i[11]) * (PESO - P5) / (P6 - P5)

    elif PESO > P6:
        if i[0] < 0.9:
            cte = -0.00043214
        elif 0.9 <= i[0] >= 1:
            cte = -0.000693476
        elif 1 < i[0]:
            cte = -0.001720016

        enve = i[1] + cte * (PESO - P6)
        alt_enve.append(enve)
        mach_enve.append(i[0])

# Los siguientes condicionales añaden el último punto a las listas
if P2 < PESO < P3:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)       
elif P3 < PESO < P4:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)
elif P4 < PESO < P5:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)
elif P5 < PESO < P6:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)

F.close() # Sentencia para cerrar el archivo de lectura "Envolventes"

# Una vez tenemos los puntos para cierto peso, realizamos una aproximación
# discreta por mínimos cuadrados
# Definimos donde empieza y termina cada tramo

# Tramo 1: (0.6 Mach,punto de altitud maxima)
c_1 = 0  # Comenzamos a comparar en 0

while alt_enve[c_1] <= alt_enve[c_1 + 1]:
    c_1 = c_1 + 1

fin_1 = c_1  # Primer punto de máxima altitud
TRAMO1_MACH = []
TRAMO1_ALT = []
for i in range(fin_1 + 1):
    TRAMO1_MACH.append(mach_enve[i])
    TRAMO1_ALT.append(alt_enve[i])

# Tramo 2: (punto de altitud maxima, punto de altitud minima)
c_2 = 0  # Volvemos a poner el contador a 0

while alt_enve[c_2 + fin_1] >= alt_enve[c_2 + fin_1 + 1]:
    c_2 = c_2 + 1

fin_2 = c_2  # Primer punto de mínima altitud
TRAMO2_MACH = []
TRAMO2_ALT = []
for i in range(fin_2 + 1):
    TRAMO2_MACH.append(mach_enve[i + fin_1])
    TRAMO2_ALT.append(alt_enve[i + fin_1])

# Tramo 3: (punto de altitud minima, mach final interpolado)
in_3 = fin_2 + fin_1
TRAMO3_MACH = []
TRAMO3_ALT = []
for i in range(size(mach_enve) - in_3):
    TRAMO3_MACH.append(mach_enve[i + in_3])
    TRAMO3_ALT.append(alt_enve[i + in_3])

TRAMOS_MACH = [TRAMO1_MACH, TRAMO2_MACH, TRAMO3_MACH]
TRAMOS_ALT = [TRAMO1_ALT, TRAMO2_ALT, TRAMO3_ALT]


def aprox_pol(x, y, k):
    """
    Realiza una aproximación por minimos cuadrados y devuelve el vector
    de coeficientes de un polinomio del grado requerido. Donde:
    x es el array de mach
    y es el array de altitudes
    k es el grado del polinomio
    Así pues obtenemos un polinomio:
        p(x) = c0 + c1·x + c2·x**2 + c3·x**3 + c4·x**4 + ···
    La función devuelve el vector de coeficientes c (array).
    """
    # Entramos con listas y las convertimos en arrays
    x = array(x)
    y = array(y)

    n = size(x)
    M = zeros((n, k + 1))

    # Calculamos la matriz M
    for i in range(n):
        for j in range(k + 1):
            M[i, j] = x[i]**j

    Mt = transpose(M)  # Matriz traspuesta de M
    A = matmul(Mt, M)
    L = cholesky(A)  # Factorización de cholesky
    b = matmul(Mt, y)
    L_inv = inv(L)  # Invertimos la matriz
    y_sol = matmul(L_inv, b)
    Lt = transpose(L)
    Lt_inv = inv(Lt)
    c = matmul(Lt_inv, y_sol)  # Matriz de coeficientes del polinomio

    return c


def coef_det(y_1, y_2, y_m):
    """
    Cálculo del coeficiente de determinación (R**2) que resulta de utilizar el
    polinomio aproximado.
    y_m = lista con valor medio de la altitud en el tramo, mismo tamaño que
          y_1 y que y_2
    y_1 = altitudes obtenidas interpolando (reales)
    y_2 = altitudes obtenidas con el polinomio (aproximadas con polinomio)
    """
    y_1 = array(y_1)
    y_2 = array(y_2)
    y_m = array(y_m)
    r_1 = y_2 - y_1  # Error residual
    r_2 = y_2 - y_m  # Variación de regresion
    R_2 = (dot(r_2, r_2))/(dot(r_1, r_1) + dot(r_2, r_2))  # Coef. determ.

    return R_2


grad = [5, 4, 6]  # Órdenes de las aproximaciones polinómicas
C = []  # Lista con 3 sublistas que son los coeficientes de los polinomios

for i in range(size(grad)):
    k = grad[i]
    P_MACH = TRAMOS_MACH[i]
    P_ALT = TRAMOS_ALT[i]

    coef = aprox_pol(P_MACH, P_ALT, k)  # Este resultado es un array
    if k == 4:
        print('Tramo', i + 1, 'p(x) =', coef[0], '+', coef[1], '·x +', coef[2],
              '·x**2 +', coef[3], '·x**3 +', coef[4], '·x**4')
    elif k == 5:
        print('Tramo', i + 1, 'p(x) =', coef[0], ' + ', coef[1], '·x +',
              coef[2], '·x**2 +', coef[3], '·x**3 +', coef[4], '·x**4 +',
              coef[5], '·x**5')
    elif k == 6:
        print('Tramo', i + 1, 'p(x) =', coef[0], ' +', coef[1], '·x +',
              coef[2], '·x**2 +', coef[3], '·x**3 +', coef[4], '·x**4 +',
              coef[5], '·x**5 +', coef[6], '·x**6')
    else:
        print('No está considerado ese grado de polinomio')

    C.append(coef)

# Cálculo de las altitudes medias para cada tramo

TRAMOS_ALT_M = []  # Lista que contiene las sublistas de altitudes medias
for i in range(3):
    ALT_M = 0
    tr_obj = TRAMOS_ALT[i]
    m = size(tr_obj)
    TRAMOS_ALT_M_i = []
    for n in range(m):
        ALT_M = ALT_M + tr_obj[n]
    ALT_M = ALT_M/size(tr_obj)  # Valor medio de la altitud en tramo i

    y_m = []
    for i in range(m):
        y_m.append(ALT_M)
    TRAMOS_ALT_M.append(y_m)

# Cálculo de las altitudes según el polinomio para cada valor de mach

TRAMOS_ALT_P = []
for i in range(3):
    co = C[i]  # Vector de coeficientes del polinomio
    s = size(co)
    mach_obj = TRAMOS_MACH[i]  # Lista de mach del tramo
    TRAMO_ALT_i = []
    for n in range(size(mach_obj)):
        ter = []  # Lista de terminos (1 x x**2 ···)
        for d in range(s):
            ap = (mach_obj[n])**d
            ter.append(ap)
        alt_p = dot(co, ter)
        TRAMO_ALT_i.append(alt_p)
    TRAMOS_ALT_P.append(TRAMO_ALT_i)

# Coeficientes de determinación de cada polinomio

R_cuad = []
for i in range(3):
    R_cuad.append(coef_det(TRAMOS_ALT[i], TRAMOS_ALT_P[i], TRAMOS_ALT_M[i]))

print('\nLos coeficientes de determinación son para cada tramo:')
print('\nTramo 1:\t', R_cuad[0])
print('\nTramo 2:\t', R_cuad[1])
print('\nTramo 3:\t', R_cuad[2])
