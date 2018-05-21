# -*- coding: utf-8 -*-
"""
@author: Team REOS

Módulo que contiene las integraciones del movimiento
"""

from numpy import inf, dot
from numpy.linalg import norm

from mecanica import aceleracion, resistencia, peso
from gravedad import RT

DT = .05


def step(mas, tie, pos, vel, gasto, isp, vloss=0, masa_minima=0, step_size=DT):
    '''Paso de integración.

    mas : float
        Masa.

    tie : float
        Tiempo.

    pos : array (3 componentes)
        Posición.

    vel : array (3 componentes)
        Velocidad.

    step_size : float
        Salto temporal. Por defecto es step_size=0.01.
    '''
    dtl = step_size
    masa = mas - gasto * dtl
    if masa < masa_minima:
        masa = masa_minima
        dtl = (mas - masa_minima) / gasto
    tiempo = tie + dtl
    acc = aceleracion(pos, vel, (mas + masa) / 2, gasto, isp)
    posicion = pos + vel * dtl + 0.5 * acc * dtl**2
    velocidad = vel + acc * dtl
    # Pérdida de velocidad
    pos_media = (pos + posicion) / 2
    vel_media = (vel + velocidad) / 2
    mas_media = (mas + masa) / 2
    loss_aero = dtl * norm(resistencia(pos_media, vel_media)) / mas_media
    loss_grav = dtl * (dot(peso(pos_media, mas_media), vel_media)
                       / (norm(vel_media) * mas_media))
    vloss = vloss + loss_aero + loss_grav

    return masa, tiempo, posicion, velocidad, vloss


def etapa(masa_etapa, masa_total, gasto, isp, posicion_inicial,
          velocidad_inicial, tiempo_inicial=0, vloss=0,
          step_size=DT, altura_maxima=inf, imprimir=False):
    '''Ejecuta todos los pasos de integración de una etapa.

    masa_inicial : float
        Masa inicial.

    pos : array (3 componentes)
        Posición inicial.

    vel : array (3 componentes)
        Velocidad inicial.

    step_size : float
        Salto temporal. Por defecto es step_size=0.01.

    altura_maxima : float
        Altura máxima. Por defecto es 500000.

    imprimir : string
        Nombre del archivo de escritura. Si False, no se escribe.
    '''
    mase = masa_etapa
    masa = masa_total
    resto = masa - mase  # Como si fuera la carga de pago
    pos = posicion_inicial
    vel = velocidad_inicial

    consumido = mase <= 0
    altur = norm(pos) - RT
    tiempo = tiempo_inicial

    while (altur < altura_maxima
           and dot(pos, vel) >= 0
           and not consumido):
        # Condiciones de parada:
        # 1) que se haya superado la altura máxima
        # 2) que el misil esté cayendo
        # 3) que se haya consumido todo el combustible de la etapa
        masa, tiempo, pos, vel, vloss = step(masa, tiempo, pos, vel, gasto,
                                             isp, vloss=vloss,
                                             masa_minima=resto,
                                             step_size=step_size)
        altur = norm(pos) - RT
        mase = masa - resto
        consumido = mase <= 0

        if imprimir:
            imprimir.write('\n' + format(tiempo, '.3f')
                           + '\t' + format(altur, '.0f')
                           + '\t' + format(norm(vel), '.0f')
                           + '\t' + format(masa, '.1f'))

    return masa, tiempo, pos, vel, vloss


def vuelo_libre(masa, posicion_inicial, velocidad_inicial, tiempo_de_vuelo=inf,
                tiempo_inicial=0, vloss=0, step_size=DT, altura_maxima=inf,
                imprimir=False):
    '''Integración del vuelo sin empuje.
    '''
    t_vuelo = 0
    tiempo = tiempo_inicial
    pos = posicion_inicial
    vel = velocidad_inicial
    altur = norm(pos) - RT
    encendido = False

    while (altur < altura_maxima
           and t_vuelo <= tiempo_de_vuelo
           and dot(pos, vel) >= 0):
        # Condiciones de parada:
        # 1) que se haya superado la altura máxima
        # 2) que se haya superado el tiempo de vuelo libre
        # 3) que el misil esté cayendo
        if tiempo_de_vuelo < step_size + t_vuelo:
            step_size = tiempo_de_vuelo - t_vuelo
            encendido = True
        masa, t_vuelo, pos, vel, vloss = step(masa, t_vuelo, pos, vel, 0, 0,
                                              vloss=vloss, step_size=step_size)
        altur = norm(pos) - RT
        tiempo = t_vuelo + tiempo_inicial
        if imprimir:
            imprimir.write('\n' + format(tiempo, '.3f')
                           + '\t' + format(altur, '.0f')
                           + '\t' + format(norm(vel), '.0f')
                           + '\t' + format(masa, '.1f'))
        if encendido:
            break

    return masa, tiempo, pos, vel, vloss


def lanzamiento(masas, estructuras, gastos, isps, posicion_inicial,
                velocidad_inicial, retardos, tiempo=0, step_size=DT,
                altura_maxima=inf, imprimir=False):
    '''Calcula el lanzamiento del cohete multietapa.
    '''
    # Condiciones iniciales
    mas = sum(masas)
    pos = posicion_inicial
    vel = velocidad_inicial
    tie = tiempo
    vloss = 0
    altur = norm(pos) - RT
    archivo = False
    if imprimir:
        archivo = open(imprimir, 'w')
        archivo.write('Tiempo (s)'
                      + '\t' + 'Altura (m)'
                      + '\t' + 'Velocidad (m/s)'
                      + '\t' + 'Masa (kg)'
                      + '\n' + format(tie, '.3f')
                      + '\t' + format(altur, '.0f')
                      + '\t' + format(norm(vel), '.0f')
                      + '\t' + format(mas, '.1f'))

    # Retardo de lanzamiento
    if retardos[0] != 0:
        mas, tie, pos, vel, vloss = vuelo_libre(mas, pos, vel,
                                                tiempo_de_vuelo=retardos[0],
                                                tiempo_inicial=tie,
                                                vloss=vloss,
                                                step_size=step_size,
                                                altura_maxima=altura_maxima,
                                                imprimir=archivo)
        altur = norm(pos) - RT
        if altur >= altura_maxima:
            if imprimir:
                archivo.close()
            return mas, tie, pos, vel, vloss
    # Primera etapa
    mas, tie, pos, vel, vloss = etapa(masas[0] * (1 - estructuras[0]), mas,
                                      gastos[0], isps[0], pos, vel,
                                      tiempo_inicial=tie, vloss=vloss,
                                      step_size=step_size,
                                      altura_maxima=altura_maxima,
                                      imprimir=archivo)
    altur = norm(pos) - RT
    if altur >= altura_maxima:
        if imprimir:
            archivo.close()
        return mas, tie, pos, vel, vloss
    mas = mas - masas[0] * estructuras[0]
    # Retardo de etapa
    if retardos[1] != 0:
        mas, tie, pos, vel, vloss = vuelo_libre(mas, pos, vel,
                                                tiempo_de_vuelo=retardos[1],
                                                tiempo_inicial=tie,
                                                vloss=vloss,
                                                step_size=step_size,
                                                altura_maxima=altura_maxima,
                                                imprimir=archivo)
        altur = norm(pos) - RT
        if altur >= altura_maxima:
            if imprimir:
                archivo.close()
            return mas, tie, pos, vel, vloss
    # Segunda etapa
    mas, tie, pos, vel, vloss = etapa(masas[1] * (1 - estructuras[1]), mas,
                                      gastos[1], isps[1], pos, vel,
                                      tiempo_inicial=tie, vloss=vloss,
                                      step_size=step_size,
                                      altura_maxima=altura_maxima,
                                      imprimir=archivo)
    altur = norm(pos) - RT
    if altur >= altura_maxima:
        if imprimir:
            archivo.close()
        return mas, tie, pos, vel, vloss
    mas = mas - masas[1] * estructuras[1]
    # Vuelo libre
    mas, tie, pos, vel, vloss = vuelo_libre(mas, pos, vel, tiempo_inicial=tie,
                                            vloss=vloss, step_size=step_size,
                                            altura_maxima=altura_maxima,
                                            imprimir=archivo)
    altur = norm(pos) - RT
    if altur >= altura_maxima:
        if imprimir:
            archivo.close()
        return mas, tie, pos, vel, vloss

    if imprimir:
        archivo.close()
    return mas, tie, pos, vel, vloss
