# -*- coding: utf-8 -*-
"""

@author: Team REOS

Modelo atmosférico: MSISE00

Funciones de temperatura (temperature), densidad (density), presión
(pressure) y viscosidad (viscosity).  Sólo requieren una variable de
entrada: la altitud.

Los datos se obtienen en la página
https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php.

"""

from scipy.interpolate import interp1d

# Constantes atmosféricas.
R_AIR = 287  # Constante de los gases ideales (J/Kkg).
RHO_SL = 101325 / (R_AIR * 288.15)  # Densidad a nivel del mar (kg/m3).
GAMMA = 1.4  # Coeficiente de dilatación adiabática.
BETA_VISC = .000001458  # Viscosidad de referencia (Pa s/K.5).
S_VISC = 110.4  # Temperatura de referencia para la viscosidad (K).

# Se importan los archivos con la temperatura y la densidad

T = open('temper', 'r')
HEIGHT = []
TEMP = []
for line in T:
    s = line.strip()
    s = s.split()
    HEIGHT.append(float(s[0]) * 1000)  # Altitude (m).
    TEMP.append(float(s[1]))  # Temperature (K).
T.close()

D = open('dens', 'r')
DENS = []
for line in D:
    s = line.strip()
    s = s.split()
    DENS.append(float(s[1]))  # Density (kg/m3).
D.close()


def temperature(alt):
    '''Cálculo de la temperatura en función de la altura dada por el
    modelo MSISE00.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 799900 metros.

    La variable de salida es un float con la temperatura (K).
    '''
    tem = interp1d(HEIGHT, TEMP)
    return float(tem(alt))


def density(alt):
    '''Cálculo de la densidad en función de la altura dada por el modelo
    MSISE00.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 799900 metros.

    La variable de salida es un float con la densidad (kg/m3).
    '''
    den = interp1d(HEIGHT, DENS)
    return float(den(alt))


def pressure(alt):
    '''Cálculo de la presión en función de la altura dada por el modelo
    MSISE00.  Se implementa la ley de los gases ideales.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 799900 metros.

    La variable de salida es un float con la presión (Pa).
    '''
    return density(alt) * R_AIR * temperature(alt)


def viscosity(alt):
    '''Cálculo de la viscosidad en función de la altura dada por el
    modelo MSISE00.  Se implementa la ley de Sutherland.

    La variable de entrada alt es la altitud (m).  Debe ser menor o
    igual que 799900 metros.

    La variable de salida es un float con la presión (Pa).
    '''
    tem = temperature(alt)
    return BETA_VISC * tem**(3 / 2) / (tem + S_VISC)
