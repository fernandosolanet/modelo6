# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:53:21 2018

@author: Team REOS

Módulo para sacar la envolvente de vuelo del avión F4 y la región optima de
lanzamiento en función de la energía mecánica.

Como input tiene el PESO del avión cargado con el REOS y la entrada manual en
la consola de la eficiencia para la energía mecánica óptima.

Se realiza una interpolación con las envolventes conocidas de unos pesos
determinados y se realiza una aproximación polinómica de alto orden para
recrear las ecuaciones de los límites de la envolvente para el PESO de entrada.
A continuación se calcula el punto de energía mecánica óptimo y una región
óptima de lanzamiento.
"""

from numpy import transpose, zeros, matmul, array, size, dot
from numpy.linalg import cholesky, inv
from numpy import linspace
import matplotlib.pyplot as plt
from Modulo_Emecanica import opt_em, region_opt, alt_min
# from inputs_iniciales import MASS

# Definción de variables. Todos los pesos en libras
# PESO = MASS # Para cuando importemos directamente del otro módulo
PESO = 39000  # Input de entrada. Peso del avión cargado con el cohete REOS
P1 = 42777.0
P2 = 43035.0
P3 = 45472.0
P4 = 46279.0
P5 = 46537.0
P6 = 48974.0

alt_enve = []  # Lista valores de altitud de la curva límite de la envolvente
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

        elif 0.9 <= i[0] <= 1:
            cte = -0.000322528

        elif 1 < i[0]:
            cte = -0.0005

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
        elif 0.9 <= i[0] <= 1:
            cte = -0.000693476
        elif 1 < i[0]:
            cte = -0.001720016

        if i[0] <= 1.5:
            enve = i[6] + cte * (PESO - P6)
            alt_enve.append(enve)
            mach_enve.append(i[0])

    elif PESO == 42777.0:
        alt_enve.append(i[1])
        mach_enve.append(i[0])

    elif PESO == 43035.0:
        alt_enve.append(i[2])
        mach_enve.append(i[0])

    elif PESO == 45472.0:
        alt_enve.append(i[3])
        mach_enve.append(i[0])

    elif PESO == 46279.0:
        alt_enve.append(i[4])
        mach_enve.append(i[0])

    elif PESO == 46537.0:
        alt_enve.append(i[5])
        mach_enve.append(i[0])

    elif PESO == 48974.0:
        alt_enve.append(i[6])
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

F.close()  # Sentencia para cerrar el archivo de levtura "Envolventes"

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

while alt_enve[c_2 + fin_1] >= alt_enve[c_2 + fin_1 + 1] or c_2 < 4:
    c_2 = c_2 + 1
    if c_2 == 8:  # Condición de salida de bucle en el caso de que no haya
        break     # mínimo y haya punto de inflexión

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
    y_m = lista del valor medio de la altitud en el tramo, mismo tamaño que
          y_1 y que y_2
    y_1 = altitudes obtenidas interpolando (reales)
    y_2 = altitudes obtenidas con el polinomio (aproximadas)
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
        print('\nTramo', i + 1, 'p(x) = {0:.2f} + {1:.2f}·x +'.format(coef[0],
                                                                      coef[1]),
              '{0:.2f}·x**2 + {1:.2f}·x**3 + {2:.2f}·x**4'.format(coef[2],
                                                                  coef[3],
                                                                  coef[4]))
    elif k == 5:
        print('\nTramo', i + 1, 'p(x) = {0:.2f} + {1:.2f}·x +'.format(coef[0],
                                                                      coef[1]),
              '{0:.2f}·x**2 + {1:.2f}·x**3 + {2:.2f}·x**4'.format(coef[2],
                                                                  coef[3],
                                                                  coef[4]),
              '{0:.2f}·x**5'.format(coef[5]))
    elif k == 6:
        print('\nTramo', i + 1, 'p(x) = {0:.2f} + {1:.2f}·x +'.format(coef[0],
                                                                      coef[1]),
              '{0:.2f}·x**2 + {1:.2f}·x**3 + {2:.2f}·x**4'.format(coef[2],
                                                                  coef[3],
                                                                  coef[4]),
              '{0:.2f}·x**5 + {1:.2f}·x**6'.format(coef[5], coef[6]))
    else:
        print('No está considerado ese grado de polinomio')

    C.append(coef)

C_1 = C[0]  # Coeficientes Tramo 1 (Grado 5)
C_2 = C[1]  # Coeficientes Tramo 2 (Grado 4)
C_3 = C[2]  # Coeficientes Tramo 3 (Grado 6)
# Calculo las altitudes medias para cada tramo

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

# Calculamos el punto óptimo

P_opt = opt_em(TRAMO3_MACH, C_3, PESO)

# Calculamos la región óptima para una eficiencia de nu = 0.85

nu = float(input('Introduzca el factor de rendimiento para EM: '))
region = region_opt(P_opt, nu, PESO, C_3, TRAMO3_MACH[-1])

x_1 = linspace(0.6, TRAMO1_MACH[-1], 50)
x_2 = linspace(TRAMO1_MACH[-1], TRAMO2_MACH[-1], 25)
x_3 = linspace(TRAMO2_MACH[-1], TRAMO3_MACH[-1], 100)

EQ_T1 = (C_1[0] + C_1[1]*x_1 + C_1[2]*x_1**2 + C_1[3]*x_1**3 + C_1[4]*x_1**4 +
         C_1[5]*x_1**5)
EQ_T2 = (C_2[0] + C_2[1]*x_2 + C_2[2]*x_2**2 + C_2[3]*x_2**3 + C_2[4]*x_2**4)
EQ_T3 = (C_3[0] + C_3[1]*x_3 + C_3[2]*x_3**2 + C_3[3]*x_3**3 + C_3[4]*x_3**4 +
         C_3[5]*x_3**5 + C_3[6]*x_3**6)

plt.plot(x_1, EQ_T1, 'k', label='Tramo 1')
plt.plot(x_2, EQ_T2, 'k', label='Tramo 2')
plt.plot(x_3, EQ_T3, 'k', label='Tramo 3')
plt.plot(P_opt[0], P_opt[2], 'ro', label='Punto de eficiencia máxima')
plt.plot(region[0], region[1], 'b.-', label='Región óptima de lanzamiento')

plt.xlabel('Mach')
plt.ylabel('Altitud ft x1000')
plt.title('Envolvente F-4 peso: {0}'.format(PESO))

plt.show()

print('La energía mecánica máxima es: {0:.2f} J'.format(P_opt[1]))
print('Para mach M = {0:.3f} y altitud {1:.3f} m'.format(P_opt[0],
                                                         P_opt[2]*1000*0.3048))
print('\nEn una parte de la región no será posible volar pues estará fuera de',
      '\nla envolvente, ya que falta el tramo de abajo de ésta, pero no es de'
      '\nnuestro inetrés pues intentaremos volar al límite de la envolvente')

M_v = float(input('Introduzca el mach de vuelo: '))
EM_opt = nu*P_opt[1]
h_min = alt_min(EM_opt, M_v, PESO)
print('Para un vuelo a '
      'M = {0} se debe volvar al menos a {1:.2f} ft'.format(M_v, h_min))
h_met = h_min*0.3048
print('En metros: {0:.2f} m'.format(h_met))
