# -*- coding: utf-8 -*-
"""

@author: TeamREOS

Funciones que nos dan los coeficientes aerodinámicos del avión F4.
Coeficientes aerodinámicos de sustentación y de resistencia y todo aquello
relacionado con los mismos. Sin empuje.

"""

from math import radians, pi, cos


#-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL VEHÍCULO-------------------

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
    '''CD0 (coeficiente de resistencia parásita).  El estudio de este
    coeficiente se da en dos partes: una incompresible y otra compresible.  Es
    función del número de Mach y de distintos coeficientes referidos a las
    partes del avión.
    '''
    #Coeficientes de las distintas partes del avión, son coeficientes
    # experimentales que se introducen en la fórmula para calcular CD0.
    x_ala = .003
    x_fuselaje = .0024
    x_gondolas = .006
    x_med = (x_fuselaje + x_gondolas) / 2
    x_cola = .0025
    #Áreas totales de las partes del avión mojdadas por la corriente de aire y
    # que, por tanto, influyen en la resistencia.
    a_ala = 2 * S_W
    a_fuselaje = .75 * pi * LF * BF
    a_gondolas = pi * LM * DM * NE
    a_cola = 2 * (S_H + S_V)
    incremento = 1.1 # Incremento por interferencias e imperfecciones (10%)
    #No se tiene en cuenta el efecto de la compresibilidad; pero su cálculo es
    # necesario para el cálculo del coeficiente de resistencia parásita con
    # compresibilidad.
    cd0_inc = ((x_ala * a_ala + x_med * a_fuselaje + x_gondolas * a_gondolas
                + x_cola * a_cola) * incremento) / S_W
    n_comp = 3 / (1 + (1 / AR))
    # Coeficiente de la fórmula con compresibilidad.
    #Cálculo de CD0 con efectos de compresibilidad.
    #Se emplea Md98, el 98% del Mach de divergencia ya que la fórmula tiende
    # asintóticamente a infinito en rangos cercanos al Mach de divergencia.
    if mach < M_D098 and mach < M_C:
        cd0_compresible = cd0_inc / ((1 - (mach / M_D)**2)**.25)
        cd_incremento = 0
        return cd0_compresible + cd_incremento
    #Para valores superiores al Mach crítico aparece un incremento por
    # resistencia de onda.
    if mach < M_D098 and mach >= M_C:
        cd0_compresible = cd0_inc / ((1 - (mach / M_D)**2)**.25)
        cd_incremento = (K / 1000) * ((10 * (mach - M_C)) / ((1 / cos(FLECHA))
                                                             - M_C))**n_comp
        return cd0_compresible + cd_incremento
    if 1.2 >= mach > M_D098:
        return (-379.32512053 * mach**5 + 1994.1499524 * mach**4
                - 4177.4704011 * mach**3 + 4358.3944768 * mach**2
                - 2264.4097020 * mach + 468.71342687)
    if mach > 1.2:
        return .031

def k(mach):
    '''Coeficiente de resistencia inducida que multiplica al coeficiente
    de sustentación.
    '''
    fos = .005 * (1 + 1.5 * (ESTRECHAMIENTO - .6)**2)
    #Este valor es una función lambda que aparece dentro del factor de Oswald.
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

def resistencia(vel, dens, c_d):
    '''Fuerza aerodinámica de resistencia total.
    '''
    return .5 * dens * S_W * c_d * vel**2

def sustentacion(vel, dens, c_l):
    '''Fuerza aerodinámica de sustentación total.
    '''
    return .5 * dens * S_W * c_l * vel**2
