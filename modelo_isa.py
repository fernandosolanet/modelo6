# -*- coding: utf-8 -*-
"""

@author: Team REOS

Modelo atmosférico: ATMÓSFERA ESTÁNDAR INTERNACIONAL (ISA).

Funciones de gradiente de temperatura (alf_isa), temperatura
(temperature), densidad (density), presión (pressure) y viscosidad
(viscosity).  Sólo requieren una variable de entrada: la altitud.

"""

from math import exp

from gravedad import MU, RT

# Constantes atmosféricas.

R_AIR = 287  # Constante de los gases ideales (J/Kkg)
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2)..
RHO_SL = 101325 / (R_AIR * 288.15)  # Densidad a nivel del mar (kg/m3).
GAMMA = 1.4  # Coeficiente de dilatación adiabática.
BETA_VISC = .000001458  # Viscosidad de referencia (Pa s/K.5).
S_VISC = 110.4  # Temperatura de referencia para la viscosidad (K).

# Estas son las alturas estipuladas según la normativa.

H_ISA1 = 11000
H_ISA2 = 20000
H_ISA3 = 32000
H_ISA4 = 47000
H_ISA5 = 51000
H_ISA6 = 71000
H_ISA7 = 84852

# Ahora se programan las variables termodinámicas, en función de la altura,
# y se relacionarán con los valores de T y alfa para cada altura estipulada.


def alfa_isa(alt):
    '''Parámetro alfa de la ISA.  Este parámetro proporciona el gradiente de
    temperatura en K/m.

    alt: altitud (m)
    '''
    if alt < H_ISA1:
        alf = -.0065
    elif alt < H_ISA2:
        alf = 0
    elif alt < H_ISA3:
        alf = .001
    elif alt < H_ISA4:
        alf = .0028
    elif alt < H_ISA5:
        alf = 0
    elif alt < H_ISA6:
        alf = -.0028
    elif alt < H_ISA7:
        alf = -.002
    else:
        alf = 0
    return alf


def temperature(alt):
    '''Cálculo de la temperatura en función de la altura dada por el
    modelo ISA.

    alt: altitud (m)
    '''
    if alt < H_ISA1:
        h_0 = 0
        t_0 = 288.15
    elif alt < H_ISA2:
        h_0 = H_ISA1
        t_0 = 216.65
    elif alt < H_ISA3:
        h_0 = H_ISA2
        t_0 = 216.65
    elif alt < H_ISA4:
        h_0 = H_ISA3
        t_0 = 228.65
    elif alt < H_ISA5:
        h_0 = H_ISA4
        t_0 = 270.65
    elif alt < H_ISA6:
        h_0 = H_ISA5
        t_0 = 270.65
    elif alt < H_ISA7:
        h_0 = H_ISA6
        t_0 = 214.65
    else:
        h_0 = H_ISA7
        t_0 = 214.65 - .002 * (H_ISA7 - H_ISA6)
    return t_0 + alfa_isa(alt) * (alt - h_0)


def density(alt):
    '''Cálculo de la densidad en función de la altura dada por el modelo
    ISA en kg/m3.  Se implementa el principio de Pascal.

    alt: altitud (m)
    '''
    t_isa = temperature(alt)

    def density_alfa_0(altit, rho_0, alt_0):
        '''Fórmula de la densidad cuando alfa es 0.
        '''
        return rho_0 * exp(-GRAV * (altit - alt_0) / (R_AIR * t_isa))

    def density_alfa_no0(rho_0, temp_0, temp_m, alfa_den):
        '''Fórmula de la densidad cuando alfa no es 0.
        '''
        return rho_0 * (temp_m / temp_0)**(-GRAV / (R_AIR * alfa_den) - 1)
    r_1 = density_alfa_no0(RHO_SL, temperature(0), temperature(H_ISA1),
                           alfa_isa(0))
    r_2 = density_alfa_0(H_ISA2, r_1, H_ISA1)
    r_3 = density_alfa_no0(r_2, temperature(H_ISA2), temperature(H_ISA3),
                           alfa_isa(H_ISA2))
    r_4 = density_alfa_no0(r_3, temperature(H_ISA3), temperature(H_ISA4),
                           alfa_isa(H_ISA3))
    r_5 = density_alfa_0(H_ISA5, r_4, H_ISA4)
    r_6 = density_alfa_no0(r_5, temperature(H_ISA5), temperature(H_ISA6),
                           alfa_isa(H_ISA5))
    r_7 = density_alfa_no0(r_6, temperature(H_ISA6), temperature(H_ISA7),
                           alfa_isa(H_ISA6))
    alf = alfa_isa(alt)
    if alt < H_ISA1:
        h_0 = 0
        t_0 = temperature(h_0)
        rho0 = RHO_SL
    if alt >= H_ISA1 and alt < H_ISA2:
        h_0 = H_ISA1
        t_0 = temperature(h_0)
        rho0 = r_1
    if alt >= H_ISA2 and alt < H_ISA3:
        h_0 = H_ISA2
        t_0 = temperature(h_0)
        rho0 = r_2
    if alt >= H_ISA3 and alt < H_ISA4:
        h_0 = H_ISA3
        t_0 = temperature(h_0)
        rho0 = r_3
    if alt >= H_ISA4 and alt < H_ISA5:
        h_0 = H_ISA4
        t_0 = temperature(h_0)
        rho0 = r_4
    if alt >= H_ISA5 and alt < H_ISA6:
        h_0 = H_ISA5
        t_0 = temperature(h_0)
        rho0 = r_5
    if alt >= H_ISA6 and alt < H_ISA7:
        h_0 = H_ISA6
        t_0 = temperature(h_0)
        rho0 = r_6
    if alt >= H_ISA7:
        h_0 = H_ISA7
        t_0 = temperature(h_0)
        rho0 = r_7
    if alf == 0:
        return density_alfa_0(alt, rho0, h_0)
    return density_alfa_no0(rho0, t_0, t_isa, alf)


def pressure(alt):
    '''Cálculo de la presión en función de la altura dada por el modelo
    ISA en Pa.  Se implementa la ley de los gases ideales.

    alt: altitud (m)
    '''

    return density(alt) * R_AIR * temperature(alt)


def viscosity(alt):
    '''Cálculo de la viscosidad en función de la altura dada por el
    modelo ISA en Pa s.  Se implementa la ley de Sutherland.

    alt: altitud (m)
    '''
    temp = temperature(alt)
    return BETA_VISC * temp**(3 / 2) / (temp + S_VISC)
