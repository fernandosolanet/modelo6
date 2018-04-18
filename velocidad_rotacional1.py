# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 14:47:11 2018

@autHor: diego.rodriguezparra
"""

from math import cos, sin
from gravedad import RT

# LATtitudes del mundo (ciudades).

# 90º Polo Norte.
# 66.57º Circulo Polar Artico.
# 50º Bruselas.
# 40º Madrid y Nueva York.
# 34º Los Angeles y Tokio.
# 27º Isla del Hierro.
# 23.43º Trópico de Cáncer.
# 19º Ciudad de México.
# 0º Ecuador Quito.
# -1º Nairobi (Kenia).
# -15.5º Brasilia.
# -16º La Paz (Bolivia).
# -23.43º Trópico de Capricornio.
# -34º Buenos Aires, Ciudad del Cabo y Sydney.
# -66.57º Circulo Polar Antartico.
# -90º Polo Sur.


def vel_rotacional(altitud, latitud, azimut):
    '''Velocidad rotacional en m/s.

    altitud :  altitud (m)
    latitud :  latitud (rad)
    azimut  :  azimut (rad)
    '''

    v_angular = .7292e-4  # Velocidad angular terretre (rad/s).

    return v_angular * (RT + altitud) * cos(latitud) * sin(azimut)  # (m/s)
