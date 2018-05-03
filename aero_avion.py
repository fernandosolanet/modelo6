# -*- coding: utf-8 -*-
"""

@author: TeamREOS

Funciones que nos dan los coeficientes aerodinámicos del avión F4.
Coeficientes aerodinámicos de sustentación y de resistencia y todo aquello
relacionado con los mismos.

"""

from math import radians, pi, cos


# ------------------CARACTERÍSTICAS GEOMÉTRICAS DEL VEHÍCULO------------------

S_W = 49.2  # Superficie alar (m2).
LF = 19.2  # Longitud del fuselaje (m).
BF = 2.87  # Envergadura del fuselaje (m).
NE = 2  # Número de motores.
DM = 0  # Diámetro de cada motor.
LM = 0  # Longitud de cada motor.
S_H = 6.39  # Área de la superficie de mando horizontal (m2).
S_V = 5.035  # Área de la superficie de mando vertical (m2).
B = 11.7  # Envergadura (m).
AR = B**2 / S_W  # Alargamiento = 2,78.
FLECHA = radians(52)  # Flecha BA.
FLECHA2 = radians(41.4)  # Flecha 1/4.
ESPESOR = .064  # Espesor relativo máximo.
K = .4  # EL perfil del F-4 es el NACA0006.4-64 por tanto la K es 0,4.
ESTRECHAMIENTO = .26  # Estrechamiento.
M_C = 1 - (.065 * (cos(FLECHA))**.6) * (100 * ESPESOR)**.637814
# Mach crítico.
M_D = 1.08 * M_C  # Mach de divergencia.
M_D098 = .98 * M_D


def cl_alfa(mach):
    '''Cálculo de la pendiente del coeficiente de sustentación, que es lineal
    respecto del ángulo de ataque.  Este coeficiente varía con respecto al Mach
    y, con el ángulo de ataque, será posible obtener el coeficiente de
    sustentación.  Como el perfil es simétrico, el término independiente de la
    recta es nulo (CL0 = 0), por lo que el coeficiente de sustentación es
    proporcional al ángulo de ataque y cl_alfa(mach) da la constante de
    proporcionalidad en función del número de Mach.
    '''
    if mach <= .8:
        return (-3.4102 * mach**5 + 1.9918 * mach**4 + 2.9597 * mach**3
                - 1.6251 * mach**2 + .4172 * mach + 2.7915)
    elif mach < 1.05:
        return 4.1216 * mach**3 - 13.25 * mach**2 + 14.343 * mach - 1.8055
    return .7534 * mach**3 - 4.2913 * mach**2 + 6.5935 * mach + .3476


def angulo_ataque(alfa_posible, mach):
    '''En función de si el ángulo de ataque obtenido es inferior o superior al
    de pérdida, la función angulo_ataque devolverá un ángulo u otro.  Cuando se
    supera el ángulo de ataque de entrada en pérdida, la función angulo_ataque
    devuelve el valor del ángulo de entrada en pérdida para así no volar en
    pérdida.
    '''
    if mach < .4:
        angulo_perdida = radians(16)
    else:
        angulo_perdida = radians(.8262 * mach**3 - .37724 * mach**2 - 6.4264
                                 * mach + 18.05)
    if alfa_posible < angulo_perdida:
        return alfa_posible
    return angulo_perdida


def cd0(mach):
    '''Esta función calcula el CD0 del avión. Las ecuaciones se han
    implementado gracias a una gráfica que representa la variación de CD0 con
    respecto al Mach para el avión  F-4 Phantom. Las instrucciones if y elif se
    refieren a los distintos tramos de la gráfica en función del número de
    Mach.
    '''
    if mach <= 0.85:

        cd_0 = 0.0277

    elif 0.85 < mach <= 1:

        cd_0 = 0.1287 * mach - 0.0817

    elif 1 < mach <= 1.05:

        cd_0 = 0.036 * mach + 0.011

    elif 1.05 < mach <= 1.20:

        cd_0 = 0.014 * mach + 0.0341

    elif 1.20 < mach <= 1.32:

        cd_0 = -0.0058 * mach + 0.0579

    elif 1.32 < mach <= 1.5:

        cd_0 = -0.0133 * mach + 0.0678

    elif 1.5 < mach <= 1.9:

        cd_0 = 0.0478

    elif 1.9 < mach <= 2:

        cd_0 = -0.019 * mach + 0.0839

    return cd_0


def k(mach):
    '''Coeficiente de resistencia inducida que multiplica al coeficiente
    de sustentación.
    '''
    fos = .005 * (1 + 1.5 * (ESTRECHAMIENTO - .6)**2)
    # Este valor es una función lambda que aparece dentro del factor de Oswald.
    e_mach = 1 / ((1 + .12 * mach**2) * (1 + (.1 * (3 * NE + 1)) / (4 + AR)
                                         + (.142 + fos * AR * (10
                                                               * ESPESOR)**.33)
                                         / (cos(FLECHA2)**2)))
    # Factor de Oswald.
    return 1 / (e_mach * pi * AR)


def cd_inducida(k_d, c_l):
    '''Coeficiente de resistencia inducida.
    '''
    return k_d * c_l**2


def cd_interferencia(mach):
    '''Coeficiente de resistencia debido a la interferencia misil-avión.
    '''
    if mach <= 0.7:
        cd_i = 0
    if 0.7 < mach < 0.955:
        cd_i = (4.9704382e3*mach**6 - 2.5004431e4*mach**5 + 5.2386095e4*mach**4
                - 5.8503607e4*mach**3 + 3.6730317e4*mach**2 - 1.2291490e4*mach
                + 1.7127795e3)
    if 0.955 < mach < 0.9655:
        cd_i = -1.6067616e1*mach**2 + 3.0869911e1*mach - 1.4799698e1
    if 0.9655 < mach < 0.99:
        cd_i = (1.963091e4*mach**4 - 7.7677873e4*mach**3 + 1.1526168e5*mach**2
                - 7.6013338e4*mach + 1.8798642e4)
    if 0.99 < mach <= 1.38:
        cd_i = 2.5938426e-2*mach**2 - 7.3453378e-2*mach + 6.6405634e-2
    if 1.38 < mach < 2:
        cd_i = 0.0002*mach**2 - 0.0008*mach + 0.0151

    return cd_i


def resistencia(vel, dens, c_d):
    '''Fuerza aerodinámica de resistencia total.
    '''
    return .5 * dens * S_W * c_d * vel**2


def sustentacion(vel, dens, c_l):
    '''Fuerza aerodinámica de sustentación total.
    '''
    return .5 * dens * S_W * c_l * vel**2
