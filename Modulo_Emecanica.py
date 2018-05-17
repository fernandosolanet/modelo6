# -*- coding: utf-8 -*-
"""
Created on Thu May  3 12:21:23 2018

@author: Team REOS

Módulo que calcula el punto óptimo de lanzamiento en función de la energía
mécanica de la que disponga el avión
"""

from modelo_msise00 import temperature
from gravedad import gravity
from numpy import array, matmul, size, linspace

# Constantes
R = 286.9
gamma = 1.4


def e_mec(P, M, h):
    '''
    Función que calcula la energía mecánica en un punto. La altura entra en
    pies y la masa en libras.
    '''
    alt = h*0.3048*1000  # Conversión de miles de pies a metros
    m = P*0.45359237  # Conversión de lb a kg
    T = temperature(alt)
    em = m*(0.5*M**2*R*gamma*T + gravity(alt)*alt)
    return em


def opt_em(T3, C_3, P):
    '''
    Inputs son el peso, la lista de Mach del tramo 3 y los coeficientes
    de ese tramo. Devuelve el punto mach con la energía mecánica máxima y
    la energía mecánica máxima.

    Devuelve un vector con 3 componentes
    - Mach donde se alcanza la EM máxima
    - El valor de la energía mecánica máxima
    - Altitud de la envolvente en ese punto
    '''
    dom = T3[-1] - T3[0]
    Min = T3[0]
    lon = int(dom/0.001)
    dM = 0.001
    EM_env = []
    mach = Min
    for i in range(lon):
        vec = []
        for n in range(size(C_3)):
            term = mach**n
            vec.append(term)
        vec = array(vec)
        C_3 = array(C_3)
        alti = matmul(vec, C_3)  # Altitud en miles de pies
        EM = e_mec(P, mach, alti)
        EM_env.append(EM)
        mach = mach + dM
    MAX_EM = max(EM_env)
    pos = EM_env.index(MAX_EM)
    Mmax = Min + pos*dM
    vec = []
    for n in range(size(C_3)):
        term = Mmax**n
        vec.append(term)
    vec = array(vec)
    C_3 = array(C_3)
    altmax = matmul(vec, C_3)
    result = [Mmax, MAX_EM, altmax]
    return result


def region_opt(EM_M, ef, P, C_3, M_max):
    '''
    Calculo de la región optima de lanzamiento. Los input son:

    EM_M: resultado de la función opt_em, es decir una lista con 2 componentes
    ef: eficiencia requerida (factor que multiplica la max Emec)
    P: Peso en libras
    C_3: array de coeficientes del polinomio del tramo 3
    M_max: Mach máximo que alcanza la envolvente
    '''
    EM_op = EM_M[1]*ef  # Energía mecánica óptima para el lanzamiento
    dM = 0.001
    EM_a = EM_M[1]  # Energía mecanica inical
    M_a = EM_M[0]  # Mach inicial
    while EM_a > EM_op:
        M_a = M_a - dM
        vec = []
        for n in range(size(C_3)):
            term = M_a**n
            vec.append(term)
        vec = array(vec)
        C_3 = array(C_3)
        alti = matmul(vec, C_3)
        EM_a = e_mec(P, M_a, alti)
    Mfin_a = M_a
    EM_d = EM_M[1]
    M_d = EM_M[0]
    while EM_d > EM_op and M_d <= M_max:
        M_d = M_d + dM
        vec = []
        for n in range(size(C_3)):
            term = M_d**n
            vec.append(term)
        vec = array(vec)
        C_3 = array(C_3)
        alti = matmul(vec, C_3)
        EM_d = e_mec(P, M_d, alti)
    Mfin_d = M_d
    num = linspace(Mfin_a, Mfin_d, 15)  # Divido en 15 el tramo óptimo
    alt_reg = []
    for v in num:  # v son los mach
        vec = []
        for n in range(size(C_3)):
            term = v**n
            vec.append(term)
        vec = array(vec)
        C_3 = array(C_3)
        alt_env = matmul(vec, C_3)
        EM_new = e_mec(P, v, alt_env)
        alt_min = alt_env
        while EM_new >= EM_op:
            alt_min = alt_min - 0.3  # Bajamos 300 pies (91.44m) en cada iter
            EM_new = e_mec(P, v, alt_min)
        alt_reg.append(alt_min)

    return [num, alt_reg]


def alt_min(EM_opt, M, P):
    '''
    Esta función nos dirá cual es la altitud mínima a la que hay que lanzar
    el cohete para cierta energía mecánica óptima.
    Los inputs son:

    - EM_opt: Energía mecánica mínima requerida para el lanzamiento
    - M: Mach de vuelo
    - P: Peso del avión
    '''
    alt = 0  # este dato es en pies, pues entra en pies para la función e_mec
    EM = e_mec(P, M, alt)  # Dato inicical
    while EM < EM_opt:
        alt = alt + 0.150  # Subimos 0.150 x1000 pies (45.72m) cada iteración
        EM = e_mec(P, M, alt)
    alt_p = alt
    result = alt_p*1000  # pasamos de miles de pies a pies
    return result
