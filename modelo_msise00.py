# -*- coding: utf-8 -*-
"""

@author: Team REOS

Modelo atmosférico: NRL-MSISE00

Funciones de temperatura (temperature), densidad (density), presión
(pressure) y viscosidad (viscosity).  Sólo requieren una variable de
entrada: la altitud.

https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php

"""


# Constantes atmosféricas.
R_AIR = 287  # Constante de los gases ideales (J/Kkg).
RHO_SL = 101325 / (R_AIR * 288.15)  # Densidad a nivel del mar (kg/m3).
GAMMA = 1.4  # Coeficiente de dilatación adiabática.
BETA_VISC = 1.458e-6  # Viscosidad de referencia (Pa s/K.5).
S_VISC = 110.4  # Temperatura de referencia para la viscosidad (K).


TEMPER = [[288.4, -.002696, -1.209e-6, 1.105e-10, -2.703e-15],  # 0-11
          [334.6, -.01953, 1.0485e-6, -2.145e-11, 1.303e-16],  # 11-20
          [225.4, -.002232, 8.624e-8, -3.209e-13],  # 20-32
          [474.95, -.00753, -6.382e-7, 3.0175e-11, -3.2195e-16],  # 32-47
          [-2828, .1816, -3.495e-6, 2.203e-11],  # 47-51
          [712.9, -.016, 1.783e-7, -7.289e-13],  # 51-71
          [1857, -.05847, 6.9745e-7, -2.8105e-12],  # 71-85
          [-4371, .1357, -1.114e-6, -2.772e-13, 2.41e-17],  # 85-105
          [17010, -.4213, 3.43e-6, -8.97e-12],  # 105-125
          [-5509, .095195, -4.916e-7, 8.732e-12],  # 125-180
          [-2638, .0488, -2.58e-7, 6.156e-13, -5.569e-19],  # 180-300
          [8.9638, .010298, -4.8383e-8, 1.1428e-13, -1.3543e-19, 6.4307e-26],
          # 300-440
          [892.6]]  # 440-800

DENSIT = [[1.2113, -1.0767e-4, 2.8143e-9],  # 0-11
          [1.144, -9.3249e-5, 2.0283e-9],  # 11-20
          [.98121, -8.5625e-5, 2.59e-9, -2.6866e-14],  # 20-32
          [.29811, -1.8891e-5, 4.0648e-10, -2.9534e-15],  # 32-47
          [.030063, -1.0345e-6, 9.0937e-6],  # 47-51
          [.028566, -1.1743e-6, 1.6299e-11, -7.6239e-17],  # 51-71
          [.0083857, -2.9077e-7, 3.3802e-12, -1.3163e-17],  # 71-85
          [.018477, -9.0294e-7, 1.7721e-11, -1.7448e-16, 8.6143e-22,
           -1.7053e-27],  # 85-105
          [5.0312e-4, -1.6648e-8, 2.0688e-13, -1.1439e-18, 2.3743e-24],
          # 105-125
          [1.4975e-5, -4.693e-10, 5.8873e-15, -3.696e-20, 1.578e-25,
           -1.4508e-31],  # 125-180
          [2.9019e-8, -4.3115e-13, 2.4264e-18, -6.1061e-24, 5.7825e-30],
          # 180-300
          [5.7268e-10, -4.077e-15, 9.7777e-21, -7.8765e-27],  # 300-440
          [1.481e-10, -1.0922e-15, 3.2335e-21, -4.7949e-27, 3.5567e-33,
           -1.0547e-39]]  # 440-800


def interval_msise00(alt):
    '''División de tramos del modelo atmosférico MSISE00

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000 (8e5) metros.
    '''
    if alt < 11000:
        i = 0
    elif alt < 20000:
        i = 1
    elif alt < 32000:
        i = 2
    elif alt < 47000:
        i = 3
    elif alt < 51000:
        i = 4
    elif alt < 71000:
        i = 5
    elif alt < 85000:
        i = 6
    elif alt < 105000:
        i = 7
    elif alt < 125000:
        i = 8
    elif alt < 180000:
        i = 9
    elif alt < 300000:
        i = 10
    elif alt < 440000:
        i = 11
    else:
        i = 12
    return i


def temperature(alt):
    '''Cálculo de la temperatura en función de la altura dada por el
    modelo MSISE00.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000(8e5) metros.

    La variable de salida es un float con la temperatura (K).
    '''
    tem = 0
    i = interval_msise00(alt)
    for j, k in enumerate(TEMPER[i]):
        tem = tem + k * alt**j
    return tem


def density(alt):
    '''Cálculo de la densidad en función de la altura dada por el modelo
    MSISE00.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000 (8e5) metros.

    La variable de salida es un float con la densidad (kg/m3).
    '''
    den = 0
    i = interval_msise00(alt)
    for j, k in enumerate(DENSIT[i]):
        den = den + k * alt**j
    return den


def pressure(alt):
    '''Cálculo de la presión en función de la altura dada por el modelo
    MSISE00.  Se implementa la ley de los gases ideales.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000 (8e5) metros.

    La variable de salida es un float con la presión (Pa).
    '''
    return density(alt) * R_AIR * temperature(alt)


def viscosity(alt):
    '''Cálculo de la viscosidad en función de la altura dada por el
    modelo MSISE00.  Se implementa la ley de Sutherland.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000 (8e5) metros.

    La variable de salida es un float con la presión (Pa).
    '''
    tem = temperature(alt)
    return BETA_VISC * tem**(3 / 2) / (tem + S_VISC)
