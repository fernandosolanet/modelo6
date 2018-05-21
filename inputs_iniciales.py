# -*- coding: utf-8 -*-
"""
@author: TeamREOS
Módulo que contiene los inputs del programa
"""

from math import radians


# Inputs de altura, posición, velocidad y masa del avión
Z0 = 12000  # Altitud inicial del misil (m)
AZ = radians(90)  # Azimut de Lanzamiento dirección
LAT = radians(0)  # Latitud inicial del avión
LON = 0  # Longitud inicial del avión
M1 = 0.8  # Número de Mach inicial
MASS = 13273 + 700  # Masa del avión cargado (kg)
N = 3.5  # Factor de carga máximo

# Caracterísiticas geométricas del misil

DIAMETRO_M = .5  # Diámetro del misil (m)
LONGITUD_CONO = .9  # Longitud del cono del misil (m)
LONGITUD_MISIL = 3  # Longitud total del misil (m)

ESPESOR_ALETA = .0065  # Espesor de la aleta (m)
CMEDIA_ALETA = .18  # Cuerda media de la aleta (m)
CRAIZ_ALETA = .24  # Cuerda raíz de la aleta (m)

SW_ALETA = .07875  # Superficie de una aleta del AIM (tomado como ref., m2)
NUM_ALETAS = 4  # Número de aletas
