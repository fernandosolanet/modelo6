# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 09:30:38 2018

@author: Team REOS

Obtención de distintas envolventes para cada 500lb (226.796kg) en un dominio de
41000 - 52500 lbs
"""
from numpy import transpose, zeros, matmul, array, size, dot
from numpy.linalg import cholesky, inv
from Modulo_envolvente_F4 import aprox_pol, coef_det

C_TRAMOS1 = []
C_TRAMOS2 = []
C_TRAMOS3 = []
for i in range(40000, 53000, 500):
    PESO = i
    P1 = 42777.0
    P2 = 43035.0
    P3 = 45472.0
    P4 = 46279.0
    P5 = 46537.0
    P6 = 48974.0

    alt_enve = []
    mach_enve = []

    F = open("Envolventes.txt", 'r')  # Apertura del fichero de lectura.

    for i in F:
        i = i.strip()
        i = i.split()
        for j, a in enumerate(i):
            i[j] = float(i[j])

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

    F.close()

    # Definimos los tramos

    # Tramo 1: (0.6 Mach,punto de altitud maxima)
    c_1 = 0

    while alt_enve[c_1] <= alt_enve[c_1 + 1]:
        c_1 = c_1 + 1

    fin_1 = c_1
    TRAMO1_MACH = []
    TRAMO1_ALT = []
    for i in range(fin_1 + 1):
        TRAMO1_MACH.append(mach_enve[i])
        TRAMO1_ALT.append(alt_enve[i])

    # Tramo 2: (punto de altitud maxima, punto de altitud minima)
    c_2 = 0

    while alt_enve[c_2 + fin_1] >= alt_enve[c_2 + fin_1 + 1] or c_2 < 4:
        c_2 = c_2 + 1
        if c_2 == 8:
            break

    fin_2 = c_2
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

    # Realizamos la aproximación discreta

    grad = [5, 4, 6]  # Órdenes de las aproximaciones polinómicas
    C = []  # Lista con 3 sublistas que son los coeficientes de los polinomios

    for i in range(size(grad)):
        k = grad[i]
        P_MACH = TRAMOS_MACH[i]
        P_ALT = TRAMOS_ALT[i]

        coef = aprox_pol(P_MACH, P_ALT, k)  # Este resultado es un array
        C.append(coef)
    C_1 = C[0]  # Array con los coeficientes del tramo 1
    C_2 = C[1]  # Array con los coeficientes del tramo 2
    C_3 = C[2]  # Array con los coeficientes del tramo 3
    C_TRAMOS1.append(C_1)
    C_TRAMOS2.append(C_2)
    C_TRAMOS3.append(C_3)

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
        R_cuad.append(coef_det(TRAMOS_ALT[i], TRAMOS_ALT_P[i],
                               TRAMOS_ALT_M[i]))
    # Escribo en el archivo de salida
    dom_1 = mach_enve[fin_1]
    dom_2 = mach_enve[in_3]
    
    OUT = open('results', 'a')
    OUT.write('Peso: {0:.1f}\n'.format(PESO))
    OUT.write('Tramo 1 (0.60, {0:.2f}): p(x)= {1:.2f}'.format(dom_1, C_1[0]))
    if C_1[1] < 0 and C_1[2] < 0:
        OUT.write(' - {0:.2f}x - {1:.2f}x**2'.format(abs(C_1[1]), abs(C_1[2])))
    elif C_1[1] < 0 and C_1[2] > 0:
        OUT.write(' - {0:.2f}x + {1:.2f}x**2'.format(abs(C_1[1]), abs(C_1[2])))
    elif C_1[1] > 0 and C_1[2] < 0:
        OUT.write(' + {0:.2f}x - {1:.2f}x**2'.format(abs(C_1[1]), abs(C_1[2])))
    else:
        OUT.write(' + {0:.2f}x + {1:.2f}x**2'.format(abs(C_1[1]), abs(C_1[2])))
    if C_1[3] < 0 and C_1[4] < 0:
        OUT.write(' - {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_1[3]),
                                                        abs(C_1[4])))
    elif C_1[3] < 0 and C_1[4] > 0:
        OUT.write(' - {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_1[3]),
                                                        abs(C_1[4])))
    elif C_1[3] > 0 and C_1[4] < 0:
        OUT.write(' + {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_1[3]),
                                                        abs(C_1[4])))
    else:
        OUT.write(' + {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_1[3]),
                                                        abs(C_1[4])))
    if C_1[5] < 0:
        OUT.write(' - {0:.2f}x**5'.format(abs(C_1[5])))
    else:
        OUT.write(' + {0:.2f}x**5'.format(abs(C_1[5])))
    OUT.write('\n\tCoef.det.: R_2= {0:.6f}\n'.format(R_cuad[0]))
    OUT.write('Tramo 2 ({0:.2f}, {1:.2f}): p(x)= '.format(dom_1, dom_2))
    if C_2[1] < 0 and C_2[2] < 0:
        OUT.write('{0:.2f} - {1:.2f}x - {2:.2f}x**2'.format(C_2[0], abs(C_2[1]),
                                                            abs(C_2[2])))
    elif C_2[1] < 0 and C_2[2] > 0:
        OUT.write('{0:.2f} - {1:.2f}x + {2:.2f}x**2'.format(C_2[0], abs(C_2[1]),
                                                            abs(C_2[2])))
    elif C_2[1] > 0 and C_2[2] < 0:
        OUT.write('{0:.2f} + {1:.2f}x - {2:.2f}x**2'.format(C_2[0], abs(C_2[1]),
                                                            abs(C_2[2])))
    else:
        OUT.write('{0:.2f} + {1:.2f}x + {2:.2f}x**2'.format(C_2[0], abs(C_2[1]),
                                                            abs(C_2[2])))
    if C_2[3] < 0 and C_2[4] < 0:
        OUT.write(' - {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_2[3]),
                                                        abs(C_2[4])))
    elif C_2[3] < 0 and C_2[4] > 0:
        OUT.write(' - {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_2[3]),
                                                        abs(C_2[4])))
    elif C_2[3] > 0 and C_2[4] < 0:
        OUT.write(' + {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_2[3]),
                                                        abs(C_2[4])))
    else:
        OUT.write(' + {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_2[3]),
                                                        abs(C_2[4])))
    OUT.write('\n\tCoef.det.: R_2= {0:.6f}\n'.format(R_cuad[1]))
    OUT.write('Tramo 3 ({0:.2f}, {1:.2f}): p(x)= {2:.2f}'.format(dom_2,
                                                                 TRAMO3_MACH[-1],
                                                                 C_3[0]))
    if C_3[1] < 0 and C_3[2] < 0:
        OUT.write(' - {0:.2f}x - {1:.2f}x**2'.format(abs(C_3[1]), abs(C_3[2])))
    elif C_3[1] < 0 and C_3[2] > 0:
        OUT.write(' - {0:.2f}x + {1:.2f}x**2'.format(abs(C_3[1]), abs(C_3[2])))
    elif C_3[1] > 0 and C_3[2] < 0:
        OUT.write(' + {0:.2f}x - {1:.2f}x**2'.format(abs(C_3[1]), abs(C_3[2])))
    else:
        OUT.write(' + {0:.2f}x + {1:.2f}x**2'.format(abs(C_3[1]), abs(C_3[2])))
    if C_3[3] < 0 and C_3[4] < 0:
        OUT.write(' - {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_3[3]),
                                                        abs(C_3[4])))
    elif C_3[3] < 0 and C_3[4] > 0:
        OUT.write(' - {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_3[3]),
                                                        abs(C_3[4])))
    elif C_3[3] > 0 and C_3[4] < 0:
        OUT.write(' + {0:.2f}x**3 - {1:.2f}x**4'.format(abs(C_3[3]),
                                                        abs(C_3[4])))
    else:
        OUT.write(' + {0:.2f}x**3 + {1:.2f}x**4'.format(abs(C_3[3]),
                                                        abs(C_3[4])))
    if C_3[5] < 0 and C_3[6] < 0:
        OUT.write(' - {0:.2f}x**5 - {1:.2f}x**6'.format(abs(C_3[5]),
                                                        abs(C_3[6])))
    elif C_3[5] < 0 and C_3[6] > 0:
        OUT.write(' - {0:.2f}x**5 + {1:.2f}x**6'.format(abs(C_3[5]),
                                                        abs(C_3[6])))
    elif C_3[5] > 0 and C_3[6] < 0:
        OUT.write(' + {0:.2f}x**5 - {1:.2f}x**6'.format(abs(C_3[5]),
                                                        abs(C_3[6])))
    else:
        OUT.write(' + {0:.2f}x**5 + {1:.2f}x**6'.format(abs(C_3[5]),
                                                        abs(C_3[6])))
    OUT.write('\n\tCoef.det.: R_2= {0:.6f}\n\n'.format(R_cuad[2]))
OUT.close()
