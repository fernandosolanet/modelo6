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


def step(mas, tie, pos, vel, gasto, isp, vloss=0, masa_minima=0, step_size=DT,
         perdidas=False):
    '''Paso de integración.

    mas : float
        Masa.

    tie : float
        Tiempo.

    pos : array (3 componentes)
        Posición.

    vel : array (3 componentes)
        Velocidad.

    gasto : float
        Gasto másico.

    isp : float
        Impulso específico.

    vloss : float
        Pérdidas de velocidad. Por defecto es vloss=0.

    masa_minima : float
        Masa mínima. Por defecto es masa_minima=0.

    step_size : float
        Salto temporal. Por defecto es step_size=DT (DT=0.05).

    perdidas : bool
        Indica si deben computarse las pérdidas de velocidad o no. Por
        defecto es perdidas=False.
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
    if perdidas:
        pos_media = (pos + posicion) / 2
        vel_media = (vel + velocidad) / 2
        mas_media = (mas + masa) / 2
        loss_aero = dtl * norm(resistencia(pos_media, vel_media)) / mas_media
        loss_grav = dtl * (dot(peso(pos_media, mas_media), vel_media)
                           / (norm(vel_media) * mas_media))
        vloss = vloss + loss_aero + loss_grav
        return masa, tiempo, posicion, velocidad, vloss
    return masa, tiempo, posicion, velocidad


def etapa(masa_etapa, masa_total, gasto, isp, posicion_inicial,
          velocidad_inicial, tiempo_inicial=0, vloss=0, step_size=DT,
          altura_maxima=inf, perdidas=False, imprimir=False):
    '''Ejecuta todos los pasos de integración de una etapa.

    masa_etapa : float
        Masa inicial de la etapa.

    masa_total : float
        Masa total inicial del lanzador.

    gasto : float
        Gasto másico.

    isp : float
        Impulso específico.

    posicion_inicial : array (3 componentes)
        Posición inicial.

    velocidad_inicial : array (3 componentes)
        Velocidad inicial.

    tiempo_inicial : float
        Tiempo inicial. Por defecto es tiempo_inicial=0.

    vloss : float
        Pérdidas de velocidad. Por defecto es vloss=0.

    step_size : float
        Salto temporal. Por defecto es step_size=DT (DT=0.05).

    altura_maxima : float
        Altura máxima. Por defecto es altura_maxima=inf.

    perdidas : bool
        Indica si deben computarse las pérdidas de velocidad o no. Por
        defecto es perdidas=False.

    imprimir : string
        Nombre del archivo de escritura. Si imprimir=False, no se
        escribe. Por defecto es imprimir=False.
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
        if perdidas:
            masa, tiempo, pos, vel, vloss = step(masa, tiempo, pos, vel, gasto,
                                                 isp, vloss=vloss,
                                                 masa_minima=resto,
                                                 step_size=step_size,
                                                 perdidas=perdidas)
        else:
            masa, tiempo, pos, vel = step(masa, tiempo, pos, vel, gasto, isp,
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
    if perdidas:
        return masa, tiempo, pos, vel, vloss

    return masa, tiempo, pos, vel


def vuelo_libre(masa, posicion_inicial, velocidad_inicial, t_de_vuelo=inf,
                tiempo_inicial=0, vloss=0, step_size=DT, altura_maxima=inf,
                perdidas=False, imprimir=False):
    '''Ejecuta todos los pasos de integración de una etapa.

    masa : float
        Masa.

    posicion_inicial : array (3 componentes)
        Posición inicial.

    velocidad_inicial : array (3 componentes)
        Velocidad inicial.

    t_de_vuelo : float

    tiempo_inicial : float
        Tiempo inicial. Por defecto es tiempo_inicial=0.

    vloss : float
        Pérdidas de velocidad. Por defecto es vloss=0.

    step_size : float
        Salto temporal. Por defecto es step_size=DT (DT=0.05).

    altura_maxima : float
        Altura máxima. Por defecto es altura_maxima=inf.

    perdidas : bool
        Indica si deben computarse las pérdidas de velocidad o no. Por
        defecto es perdidas=False.

    imprimir : string
        Nombre del archivo de escritura. Si imprimir=False, no se
        escribe. Por defecto es imprimir=False.
    '''
    t_vuelo = 0
    tiempo = tiempo_inicial
    pos = posicion_inicial
    vel = velocidad_inicial
    altur = norm(pos) - RT
    encendido = False

    while (altur < altura_maxima
           and t_vuelo <= t_de_vuelo
           and dot(pos, vel) >= 0):
        # Condiciones de parada:
        # 1) que se haya superado la altura máxima
        # 2) que se haya superado el tiempo de vuelo libre
        # 3) que el misil esté cayendo
        if t_de_vuelo < step_size + t_vuelo:
            step_size = t_de_vuelo - t_vuelo
            encendido = True
        if perdidas:
            masa, t_vuelo, pos, vel, vloss = step(masa, t_vuelo, pos, vel, 0,
                                                  0, vloss=vloss,
                                                  step_size=step_size,
                                                  perdidas=perdidas)
        else:
            masa, t_vuelo, pos, vel = step(masa, t_vuelo, pos, vel, 0, 0,
                                           step_size=step_size)
        altur = norm(pos) - RT
        tiempo = t_vuelo + tiempo_inicial
        if imprimir:
            imprimir.write('\n' + format(tiempo, '.3f')
                           + '\t' + format(altur, '.0f')
                           + '\t' + format(norm(vel), '.0f')
                           + '\t' + format(masa, '.1f'))
        if encendido:
            break

    if perdidas:
        return masa, tiempo, pos, vel, vloss

    return masa, tiempo, pos, vel


def lanzamiento(masas, estructuras, gastos, isps, posicion_inicial,
                velocidad_inicial, retardos, tiempo=0, step_size=DT,
                alt_maxima=inf, perdidas=False, imprimir=False):
    '''Ejecuta todos los pasos de integración de una etapa.

    masas : array
        Masas de las distintas etapas, incluyendo la masa de la carga de
        pago como último elemento.

    estructuras : array
        Razones estructurales de las distintas etapas.

    gastos : array
        Gastos másicos de las distintas etapas.

    isps : array
        Impulsos específicos de las distintas etapas.

    posicion_inicial : array (3 componentes)
        Posición inicial.

    velocidad_inicial : array (3 componentes)
        Velocidad inicial.

    retardos : array
        Tiempos de retardo de encendido de las distintas etapas.

    tiempo : float
        Tiempo inicial. Por defecto es tiempo=0.

    step_size : float
        Salto temporal. Por defecto es step_size=DT (DT=0.05).

    altura_maxima : float
        Altura máxima. Por defecto es altura_maxima=inf.

    perdidas : bool
        Indica si deben computarse las pérdidas de velocidad o no. Por
        defecto es perdidas=False.

    imprimir : string
        Nombre del archivo de escritura. Si imprimir=False, no se
        escribe. Por defecto es imprimir=False.
    '''
    # Condiciones iniciales
    mas = sum(masas)
    pos = posicion_inicial
    vel = velocidad_inicial
    tie = tiempo
    if perdidas:
        per = 0
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

    for i, gas in enumerate(gastos):
        # Retardo de encendido
        if retardos[i] != 0:
            if perdidas:
                mas, tie, pos, vel, per = vuelo_libre(mas, pos, vel,
                                                      t_de_vuelo=retardos[i],
                                                      tiempo_inicial=tie,
                                                      vloss=per,
                                                      step_size=step_size,
                                                      altura_maxima=alt_maxima,
                                                      perdidas=perdidas,
                                                      imprimir=archivo)
            else:
                mas, tie, pos, vel = vuelo_libre(mas, pos, vel,
                                                 t_de_vuelo=retardos[i],
                                                 tiempo_inicial=tie,
                                                 step_size=step_size,
                                                 altura_maxima=alt_maxima,
                                                 imprimir=archivo)
            altur = norm(pos) - RT
            if altur >= alt_maxima:
                if imprimir:
                    archivo.close()
                if perdidas:
                    return mas, tie, pos, vel, per
                return mas, tie, pos, vel
        # Primera etapa
        if perdidas:
            mas, tie, pos, vel, per = etapa(masas[i] * (1 - estructuras[i]),
                                            mas, gas, isps[i], pos, vel,
                                            tiempo_inicial=tie, vloss=per,
                                            step_size=step_size,
                                            altura_maxima=alt_maxima,
                                            perdidas=perdidas,
                                            imprimir=archivo)
        else:
            mas, tie, pos, vel = etapa(masas[i] * (1 - estructuras[i]), mas,
                                       gas, isps[i], pos, vel,
                                       tiempo_inicial=tie, step_size=step_size,
                                       altura_maxima=alt_maxima,
                                       imprimir=archivo)
        altur = norm(pos) - RT
        if altur >= alt_maxima:
            if imprimir:
                archivo.close()
            if perdidas:
                return mas, tie, pos, vel, per
            return mas, tie, pos, vel
        mas = mas - masas[i] * estructuras[i]
    # Vuelo libre
    if perdidas:
        mas, tie, pos, vel, per = vuelo_libre(mas, pos, vel,
                                              tiempo_inicial=tie, vloss=per,
                                              step_size=step_size,
                                              altura_maxima=alt_maxima,
                                              perdidas=perdidas,
                                              imprimir=archivo)
    else:
        mas, tie, pos, vel = vuelo_libre(mas, pos, vel, tiempo_inicial=tie,
                                         step_size=step_size,
                                         altura_maxima=alt_maxima,
                                         imprimir=archivo)

    if imprimir:
        archivo.close()
    if perdidas:
        return mas, tie, pos, vel, per
    return mas, tie, pos, vel
