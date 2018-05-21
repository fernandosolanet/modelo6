# -*- coding: utf-8 -*-
"""
Created on Thu May  3 09:16:09 2018

@author: fernando.solanet

La idea es calcular, vectorialmente, las fuerzas que actúan sobre el
misil. Con las fuerzas calculadas, se obtiene la aceleración. Con la
aceleración, se calcula la velocidad en el punto siguiente. Con la
velocidad y la aceleración, se calcula la posición en el punto
siguiente.

Falta tener en cuenta la masa estructural.
"""

from time import time
from numpy import array, pi, cross, radians, degrees, inf, exp, sqrt
from numpy.linalg import norm
from scipy.optimize import brenth

from gravedad import vel_orbital, RT
from inputs_iniciales import LAT, LON, AZ
from velocidad_rotacional1 import OMEGA_R
from avion import DT as DT_AVION
from mecanica import energia_mecanica, G0
from apoyo import angulo_vectores, inicio, cartesianas, vector_esf
from integracion import lanzamiento

TIME = time()

GASTOS = array([14, 10])  # Gasto másico (kg/s)
ISPS = array([280, 280])  # Impulso específico (s)
MASA_DE_PAGO = 143  # Masa de la carga de pago (kg)
MASAS = array([384, 173, MASA_DE_PAGO])
ESTRUCTURAS = array([0.2, 0.2])  # Razón estructural
RETARDOS = array([10, 5])
# Tiempos de retardo de encendido de las etapas (s)
DT = .05

EXPORTS = open('exports', 'r')
ALTURA = EXPORTS.readline().split()
VELOCIDAD = EXPORTS.readline().split()
FI_LIST = EXPORTS.readline().split()
PSI_LIST = EXPORTS.readline().split()
EXPORTS.close()


def condiciones_iniciales(indice, latitud, longitud, azimut):
    '''Define las condiciones iniciales (tiempo <tiempo>, posición
    <posicion> y velocidad <velocidad>) de cada punto de la trayectoria
    del avión.

    indice : int
        Índice que señala cada punto de la trayectoria del avión.

    latitud : float
        Latitud inicial del avión (rad).

    longitud : float
        Longitud inicial del avión (rad).

    azimut : float
        Azimut inicial del avión (rad).
    '''
    tiempo = indice * DT_AVION
    altura = float(ALTURA[indice])
    psi = radians(float(PSI_LIST[indice]))
    velocidad_inicial = float(VELOCIDAD[indice])
    fil = radians(float(FI_LIST[indice]))

    r_inicial = RT + altura
    s_inicial = psi * RT
    theta, phi, alpha = inicio(pi / 2 - latitud, longitud, azimut, s_inicial)
    posicion = cartesianas([r_inicial, theta, phi])
    velocidad = (vector_esf(velocidad_inicial, fil, alpha, posicion)
                 + cross(OMEGA_R, posicion))

    return tiempo, posicion, velocidad


MAXIMO = -inf
INDICE_MAXIMO = 0

print(format(sum(MASAS), '.2f') + ' kg')

for i in enumerate(ALTURA):
    t, x, v = condiciones_iniciales(i[0], LAT, LON, AZ)

    m, t, x, v, vloss = lanzamiento(MASAS, ESTRUCTURAS, GASTOS, ISPS, x, v,
                                    RETARDOS, tiempo=t, altura_maxima=7.9e5)
    energia = energia_mecanica(m, x, v)

    if MAXIMO < energia:
        MAXIMO = energia
        INDICE_MAXIMO = i[0]

    print('\r' + format(i[0] / (len(ALTURA) - 1), '.1%') + ' completado  ',
          end='\r')

T, X, V = condiciones_iniciales(INDICE_MAXIMO, LAT, LON, AZ)
V0 = norm(V)
M, TIEMPO, X, V, VLOSS = lanzamiento(MASAS, ESTRUCTURAS, GASTOS, ISPS, X, V,
                                     RETARDOS, tiempo=T, step_size=DT/5,
                                     altura_maxima=7.999e5, imprimir='14')

# Optimizar masas
i = ISPS[1] / ISPS[0]
r1 = ESTRUCTURAS[0]
r2 = ESTRUCTURAS[1]
delta = norm(V) - V0 + VLOSS
c = G0 * ISPS[0]
k = (1 + r1) * (1 + r2)**i * exp(-delta / c)


def eq(x_solve):
    '''Encontrar f2.
    '''
    if i == 1:
        return x_solve - sqrt(k * r2 / r1) + r2
    return (x_solve + r2)**(i + 1) / (r2 - x_solve * (i - 1)) - k / r1


f2 = brenth(eq, 0, 1)
f1 = i * f2 * r1 / (r2 - f2 * (i - 1))

MASA1 = sum(MASAS) * (1 - f1)
MASA2 = sum(MASAS) * f1 * (1 - f2)
MASA_UTIL = sum(MASAS) - MASA1 - MASA2

print('\n\nÁngulo de asiento: '
      + format(90 - degrees(angulo_vectores(X, V)), '.2f') + ' deg')
print('Velocidad del misil: '
      + format(norm(V) / vel_orbital(norm(X) - RT), '.3%')
      + ' de la velocidad orbital')
print('Energía mecánica: ' + format(MAXIMO, '.3e') + ' J')

print('\nMasa final: ' + format(M, '.2f') + ' kg')
print('Tiempo de lanzamiento: ' + format(TIEMPO, '.3f') + ' s')
print('Altura final: ' + format(norm(X) - RT, '.0f') + ' m')
print('Velocidad final: ' + format(norm(V), '.2f') + ' m/s')

print('\nMasa de la primera etapa: ' + format(MASA1, '.0f') + ' kg'
      + '\nMasa de la segunda etapa: ' + format(MASA2, '.0f') + ' kg'
      + '\nMasa útil: ' + format(MASA_UTIL, '.0f') + ' kg')

print('\nTiempo de ejecución: ' + format(time() - TIME, '.4f') + ' s')
