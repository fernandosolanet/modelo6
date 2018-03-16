# -*- coding: utf-8 -*-
"""
@author: Team REOS

Modelo gravitatorio newtoniano.

"""

# Constantes gravitatorias
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).


def gravity(alt):
    '''Cálculo de la gravedad en función de la altura en m/s2.  Se
    implementa la ley de gravitación universal de Newton.

    alt: altitud (m)
    '''
    radio = RT + alt
    return MU / radio**2
