# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 13:13:16 2018

@author: pedro.munoz
"""

from math import exp, radians, cos, pi, sin, degrees, atan, log10
from modeloISA import density, temperature, GAMMA, viscosity, pressure, R_AIR

#---------------------------CURVAS DE EMPUJE DEL F-18---------------------------
#Funciones de las curvas de empuje de los motores del F-18. Modelo GE F404.
#Previamente a este código, se han obtenido unas funciones de coeficientes
#relativas a las curvas de empuje de la aeronave.  Partiendo del empuje máximo
#del F-18 (el empuje máximo a nivel del mar) y con el dato de la altura, el cual
#nos dará el valor de densidad, podremos obtener el valor del empuje equivalente
#a cada altura en la que nos estemos moviendo.

TH_SL = 100000  # Empuje a nivel del mar (máximo).

def thrust(mach, den):
    '''Cálculo del empuje de la aeronave. Esta función se obtiene a partir de
    las gráficas del empuje del motor GE F404-400, cf. "Thrust Data for
    Performance Calculations" en M. Saarlas, "Aircraft Performance", p. 273.
    Se tiene en cuenta que la aeronave cuenta con dos motores GE F404-400.
    '''
    
    d_th = den / density(0)
    i = (.050618013228 + .11323534299 * d_th + 7.8263530571 * d_th**2
         - 15.012158645 * d_th**3)
    a_th = 1.1062543547 * d_th**1.276913816
    c_th = d_th * .862301392 + 1.937299323
    z_th = -.347382668*d_th + 1.71160358
    return TH_SL * (a_th + i * exp(-c_th * (mach - z_th)**2))
