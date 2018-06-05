# -*- coding: utf-8 -*-
"""
@author: Team REOS
Modelo atmosférico: NRL-MSISE00
Funciones de temperatura (temperature), densidad (density), presión
(pressure) y viscosidad (viscosity).  Sólo requieren una variable de
entrada: la altitud.
Los datos se obtienen de la página
https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
"""

from errores import ValorInadmisibleError

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
          [-5498.54, .0950362, -4.90839e-7, 8.72167e-13],  # 125-180
          [-2638, .0488, -2.58e-7, 6.156e-13, -5.569e-19],  # 180-300
          [8.9638, .010298, -4.8383e-8, 1.1428e-13, -1.3543e-19, 6.4307e-26],
          # 300-440
          [892.6]]  # 440-800

DENSIT = [[1.2113, -1.0767e-4, 2.8143e-9],  # 0-11
          [1.144, -9.3249e-5, 2.0283e-9],  # 11-20
          [.98121, -8.5625e-5, 2.59e-9, -2.6866e-14],  # 20-32
          [.53696369783, -4.3597202128e-5, 1.3582294965e-9, -1.9138461925e-14,
           1.0253010619e-19],  # 32-47
          [.030063, -1.0345e-6, 9.0937e-12],  # 47-51
          [.028566, -1.1743e-6, 1.6299e-11, -7.6239e-17],  # 51-71
          [.021858322344, -9.8512966001e-7, 1.6779533018e-11,
           -1.279121248e-16, 3.6795165237e-22],  # 71-85
          [.018477335855, -9.0294288939e-7, 1.7720644879e-11,
           -1.7447869942e-16, 8.6142822069e-22, -1.7052994766e-27],  # 85-105
          [2.1772683977e-3, -8.9767045236e-8, 1.4829610326e-12,
           -1.2267775208e-17, 5.0809300164e-23, -8.427136289e-29],  # 105-125
          [6.21790234479e-5, -2.35674472104e-9, 3.72303267358e-14,
           -3.13620341158e-19, 1.48529238085e-24, -3.74874972816e-30,
           3.938525302088e-36],  # 125-180
          [1.6244922565e-7, -3.6911722938e-12, 3.530199049e-17,
           -1.8124348491e-22, 5.2555686765e-28, -8.1472260337e-34,
           5.2685895499e-40],  # 180-300
          [3.2983479421e-9, -3.9057260867e-14, 1.8766393702e-19,
           -4.5571029771e-25, 5.5779327744e-31, -2.7476359591e-37],  # 300-440
          [3.17056098e-10, -2.79591705e-15, 1.03325853e-20, -2.04426659e-26,
           2.28016783e-32, -1.35785279e-38, 3.36960028e-45]]  # 440-800


def interval_msise00(alt):
    '''División de tramos del modelo atmosférico MSISE00.
    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 800000 (8e5) metros.
    '''
    if alt < 0:
        raise ValorInadmisibleError(dict(alt=alt), '.0f', 'positivo')
    elif alt > 8e5:
        raise ValorInadmisibleError(dict(alt=alt), '.0f', 'menor que 800000')
    i = -1
    tramos = [0, 11e3, 20e3, 32e3, 47e3, 51e3, 71e3, 85e3, 105e3, 125e3, 180e3,
              300e3, 440e3, 800e3]
    en_tramo = False
    while not en_tramo:
        i = i + 1
        en_tramo = tramos[i] <= alt and alt <= tramos[i+1]
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
