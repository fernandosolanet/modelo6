# -*- coding: utf-8 -*-
"""
@author: Team REOS

Este módulo contiene la mecánica del lanzamiento.
"""

from numpy import sqrt, cross, dot
from numpy.linalg import norm

from gravedad import gravity, MU, RT
from modelo_msise00 import temperature, pressure, GAMMA, R_AIR
from velocidad_rotacional1 import OMEGA_R
from aero_misil import cdll, SREF_MISIL

G0 = 9.81  # Constante de normalización del impulso específico (m/s2)


def empuje(gasto, impulso, vel):
    '''Calcula el empuje del lanzador.

    gasto : float
        Gasto másico.

    impulso : float
        Impulso específico normalizado.

    vel : array (3 componentes)
        Vector velocidad.
    '''
    return gasto * G0 * impulso * vel / norm(vel)


def numero_mach(pos, vel):
    '''Calcula el número de Mach.

    pos : array (3 componentes)
        Vector posición.

    vel : array (3 componentes)
        Vector velocidad.
    '''
    altur = norm(pos) - RT
    tem = temperature(altur)
    vel_sonido = sqrt(GAMMA * R_AIR * tem)
    vel_aire = cross(OMEGA_R, pos)
    vel_relativa = vel - vel_aire
    return norm(vel_relativa) / vel_sonido


def resistencia(pos, vel):
    '''Calcula la fuerza de resistencia del lanzador.

    pos : array (3 componentes)
        Vector posición.

    vel : array (3 componentes)
        Vector velocidad.
    '''
    altur = norm(pos) - RT
    mach = numero_mach(pos, vel)
    cd_misil = cdll(mach, altur)

    return -(.5 * GAMMA * pressure(altur) * mach**2 * SREF_MISIL * cd_misil
             * vel / norm(vel))


def peso(pos, mas):
    '''Calcula el peso del lanzador.

    pos : array (3 componentes)
        Vector posición.

    mas : float
        Masa.
    '''
    altur = norm(pos) - RT

    return -mas * gravity(altur) * pos / norm(pos)


def aceleracion(pos, vel, mas, gasto, isp):
    '''Calcula la aceleración total del lanzador.

    pos : array (3 componentes)
        Vector posición.

    vel : array (3 componentes)
        Vector velocidad.

    mas : float
        Masa.

    gasto : float
        Gasto másico.

    isp : float
        Impulso específico.
    '''
    emp = empuje(gasto, isp, vel)  # Empuje
    res = resistencia(pos, vel)  # Resistencia
    pes = peso(pos, mas)  # Peso

    fuerza = emp + res + pes  # Fuerza resultante total

    return fuerza / mas


def energia_mecanica(mas, pos, vel):
    '''Calcula la energía mecánica del lanzador.

    mas : float
        Masa.

    pos : array (3 componentes)
        Vector posición.

    vel : array (3 componentes)
        Vector velocidad.
    '''
    cin = .5 * mas * dot(vel, vel)
    pot = -MU * mas / norm(pos)
    mec = cin + pot

    return mec
