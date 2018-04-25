# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 16:20:30 2018

@author: ignacio.garciaguerre
"""

from math import degrees

from modelo_msise00 import density, temperature, GAMMA, viscosity, pressure, R_AIR
from aero_avion import cl_alfa, angulo_ataque, S_W
from gravedad import gravity
from modelo_empuje_pc import thrust
from Envolvente_4 import MASS, altura_techo

# -------------------CONDICIONES GRAVITATORIAS-------------------------
G = 6.673e-11  # Constante de gravitacióN universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # AceleracióN de la gravedad a nivel del mar (m/s2).
W = MASS * GRAV

for altura in range(8000, 17000, 1000):
    M = 0.6
    alturatexto = str(altura)
    F = open(alturatexto, "w")  # Fichero de escritura sin extensióN.
    F.write('MACH \tTecho \tAltura \tCD\tAngulo perdida \tANGULO DE ATAQUE\n')
    while M < 2.04:

        techo = altura_techo(M)  # valor en ft
        techo_metros = techo * 0.3048 * 1000
        rho = density(altura)
        T = temperature(altura)
        Empuje = thrust(M, rho)
        v = M * (GAMMA * R_AIR * T)**0.5

        # Ecuación eje horizontal en crucero
        CD = 2 * Empuje / (S_W * rho * v**2)

        # Ecuacion eje vertical
        CL = 2 * W / (S_W * rho * v**2)
        # ángulo de ataque
        CL_alfa = cl_alfa(M)
        alfa = CL / CL_alfa
        alfa_grados = degrees(alfa)
        alfa_perdida_actual = angulo_ataque(alfa, M)
        perdida = alfa - alfa_perdida_actual

        F.write('%.8f\t' % M)
        F.write('%.8f\t' % techo_metros)
        F.write('%.8f\t' % altura)
        F.write('%.8f\t' % CD)
        F.write('%.8f\t' % CL_alfa)
        F.write('%.8f\n' % alfa_grados)

        M = M + 0.02
        
F.close()