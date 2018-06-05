# -*- coding: utf-8 -*-
"""
@author: Team REOS

Este módulo contiene las funciones matemáticas necesarias para la
función <lanzamiento> del módulo integracion.py.
"""

from numpy import (arctan, sqrt, pi, sign, array, cos, sin, arctan2, radians,
                   cross)
from numpy.linalg import norm

import sys
sys.path.insert(0, '/path/to/modulos')
sys.path.insert(0, '/path/to/aerodinamica')
sys.path.insert(0, '/path/to/atmosfera')
sys.path.insert(0, '/path/to/empuje')

from modulos.atmosfera.gravedad import RT
from modulos.velocidad_rotacional1 import OMEGA_R
from avion import DT as DT_AVION


def esfericas(vector):
    '''Calcula las coordenadas esféricas de un vector de posición dado
    en coordenadas cartesianas.

    vector : array (3 componentes)
        vector[0] : float
            Componente x.
        vector[1] : float
            Componente y.
        vector[2] : float
            Componente z.
    '''
    # Cartesianas
    x_esf = vector[0]
    y_esf = vector[1]
    z_esf = vector[2]

    # Esféricas
    r_esf = norm(vector)

    if z_esf < 0:
        theta_esf = arctan(sqrt(x_esf**2 + y_esf**2) / z_esf) + pi
    elif z_esf == 0:
        theta_esf = pi / 2
    else:
        theta_esf = arctan(sqrt(x_esf**2 + y_esf**2) / z_esf)

    if x_esf > 0 and y_esf > 0:
        phi_esf = arctan(y_esf / x_esf)
    elif x_esf > 0 and y_esf < 0:
        phi_esf = arctan(y_esf / x_esf) + 2 * pi
    elif x_esf == 0:
        phi_esf = pi * sign(y_esf) / 2
    else:
        phi_esf = arctan(y_esf / x_esf) + pi

    return array([r_esf, theta_esf, phi_esf])


def cartesianas(vector):
    '''Calcula las coordenadas cartesianas de un vector de posición dado
    en coordenadas esféricas.

    vector : array (3 componentes)
        vector[0] : float
            Componente r. r > 0.
        vector[1] : float
            Componente theta (rad). 0 <= theta <= pi.
        vector[2] : float
            Componente phi (rad). 0 <= phi < 2*pi.
    '''
    # Cartesianas
    r_car = vector[0]
    theta_car = vector[1]
    phi_car = vector[2]

    # Esféricas
    x_car = r_car * round(sin(theta_car), 15) * round(cos(phi_car), 15)
    y_car = r_car * round(sin(theta_car), 15) * round(sin(phi_car), 15)
    z_car = r_car * round(cos(theta_car), 15)

    return array([x_car, y_car, z_car])


def inicio(theta0, phi0, alpha0, arc0):
    '''Calcula la colatitud <theta1>, la longitud <phi1> y el azimut
    <alpha1> del punto inicial de la trayectoria del lanzador.

    theta0 : float
        Colatitud inicial (rad).

    phi0 : float
        Longitud inicial (rad).

    alpha0 : float
        Azimut inicial (rad).

    arc0 : float
        Distancia recorrida.
    '''
    sin_alpha = sin(theta0) * sin(alpha0)
    sigma = arc0 / RT

    latit = arctan2(cos(theta0) * cos(sigma)
                    + sin(theta0) * sin(sigma) * cos(alpha0),
                    sqrt(sin_alpha**2
                         + (cos(theta0) * sin(sigma)
                            - sin(theta0) * cos(sigma) * cos(alpha0))**2))
    theta1 = pi / 2 - latit

    delta = arctan2(sin(alpha0) * sin(sigma),
                    sin(theta0) * cos(sigma)
                    - cos(theta0) * sin(sigma) * cos(alpha0))
    phi1 = phi0 + delta

    alpha1 = arctan2(sin_alpha,
                     -cos(theta0) * sin(sigma)
                     + sin(theta0) * cos(sigma) * cos(alpha0))

    return theta1, phi1, alpha1


def vector_esf(modulo, inc, azimut, pos):
    '''Calcula las componentes cartesianas de un vector.

    modulo : float
        Módulo del vector. modulo > 0.

    inc : float
        Ángulo sobre la horizontal local (rad). -pi/2 <= inc <= pi/2.

    azimut : float
        Azimut (rad). -pi < azimut <= pi.

    pos : array (3 componentes)
        Vector posición en coordenadas cartesianas.
    '''
    posi = esfericas(pos)
    thet = posi[1]
    phip = posi[2]
    r_vec = modulo * round(sin(inc), 15)
    t_vec = modulo * round(cos(inc), 15) * round(cos(azimut), 15)
    p_vec = modulo * round(cos(inc), 15) * round(sin(azimut), 15)

    x_1 = (r_vec * round(sin(thet), 15) * round(cos(phip), 15)
           + t_vec * round(cos(thet), 15) * round(cos(phip), 15)
           - p_vec * round(sin(phip), 15))
    x_2 = (r_vec * round(sin(thet), 15) * round(sin(phip), 15)
           + t_vec * round(cos(thet), 15) * round(sin(phip), 15)
           + p_vec * round(cos(phip), 15))
    x_3 = (r_vec * round(cos(thet), 15)
           - t_vec * round(sin(thet), 15))

    return array([x_1, x_2, x_3])


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
