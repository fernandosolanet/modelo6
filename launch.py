# -*- coding: utf-8 -*-
"""
Created on Thu May  3 09:16:09 2018

@author: fernando.solanet
"""

from time import time
from numpy import array, inf
from numpy.linalg import norm

from inputs_iniciales import LAT, LON, AZ
from gravedad import vel_orbital, RT
from mecanica import energia_mecanica
from apoyo import condiciones_iniciales, ALTURA
from integracion import lanzamiento

TIME = time()

GASTOS = array([20, 20])  # Gasto másico (kg/s)
ISPS = array([280, 280])  # Impulso específico (s)
MASA_DE_PAGO = 40  # Masa de la carga de pago (kg)
MASAS = array([635, 215, MASA_DE_PAGO])
MASA_TOTAL = float(sum(MASAS))
ESTRUCTURAS = array([0.2, 0.2])  # Razón estructural
RETARDOS = array([10, 5])
# Tiempos de retardo de encendido de las etapas (s)
DT = .05

if __name__ == '__main__':
    MAXIMO = -inf
    INDICE_MAXIMO = 0

    print(format(sum(MASAS), '.2f') + ' kg')

    for i in enumerate(ALTURA):
        t, x, v = condiciones_iniciales(i[0], LAT, LON, AZ)

        m, t, x, v = lanzamiento(MASAS, ESTRUCTURAS, GASTOS, ISPS, x, v,
                                 RETARDOS, tiempo=t, alt_maxima=7.999e5)
        energia = energia_mecanica(m, x, v)

        if MAXIMO < norm(x) - RT:
            MAXIMO = norm(x) - RT
            INDICE_MAXIMO = i[0]

        print('\r' + format(i[0] / (len(ALTURA) - 1), '.1%') + ' completado  ',
              end='\r')

    T, X, V = condiciones_iniciales(INDICE_MAXIMO, LAT, LON, AZ)
    V0 = norm(V)
    M, TIEMPO, X, V, VLOSS = lanzamiento(MASAS, ESTRUCTURAS, GASTOS, ISPS, X,
                                         V, RETARDOS, tiempo=T, step_size=DT/9,
                                         alt_maxima=7.999e5, perdidas=True,
                                         imprimir='14')

#    print('\n\nÁngulo de asiento: '
#          + format(90 - degrees(angulo_vectores(X, V)), '.2f') + ' deg')
    print('Velocidad del misil: '
          + format(norm(V) / vel_orbital(norm(X) - RT), '.3%')
          + ' de la velocidad orbital')
#    print('Energía mecánica: ' + format(energia_mecanica(M, X, V), '.3e')
#          + ' J')

#    print('\nMasa final: ' + format(M, '.2f') + ' kg')
#    print('\nTiempo de lanzamiento: ' + format(TIEMPO, '.3f') + ' s')
    print('Altura final: ' + format((norm(X) - RT) / 1000, '.3f') + ' km')

    print('\nMasa de la primera etapa: ' + format(MASAS[0], '.0f') + ' kg'
          + '\nMasa de la segunda etapa: ' + format(MASAS[1], '.0f') + ' kg'
          + '\nMasa de la tercera etapa: ' + format(MASAS[2], '.0f') + ' kg'
          + '\nMasa útil: ' + format(MASAS[3], '.0f') + ' kg')

    print('\nTiempo de ejecución: ' + format(time() - TIME, '.4f') + ' s')
