# -*- coding: utf-8 -*-
"""
@author: Team REOS
"""

from math import radians, cos, sin, degrees, pi
from modeloISA import density, temperature, GAMMA, viscosity, pressure, R_AIR
from gravedad import gravity
from modelo_empuje import thrust
from aero_avion_interferencia import cl_alfa, angulo_ataque, k, cd0
from aero_avion_interferencia import cd_inducida, S_W, CD_interferencia
from aero_avion_interferencia import resistencia, sustentacion
from aero_misil import cdll, SREF_MISIL
from velocidad_crucero_interferencia import vuelo_crucero


# --------------------------CONDICIONES GRAVITATORIAS--------------------------
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2).


# ------------------------CARACTERÍSTICAS DE LA AERONAVE-----------------------

N = 3.5  # Factor de carga máximo.
MASS = 14273  # Masa del avión cargado (kg).
W = MASS * GRAV  # Peso del avión (N).

BETA = 89
BETA = radians(BETA)
F = open('resultados', 'w')  # Fichero de escritura sin extensión.

# Las variables que se escriben a continuación corresponden al paso de
# integración del avión en cada paso
F.write('TIME (s)\tALTURA lanzamiento (m) \tDesplazamiento (m)')
F.write('\tVELOCIDAD (m/s)\tMACH\tCD_avion \tRESISTENCIA')
F.write('\tEMPUJE \tALFA (deg) \tGAMMA (deg)')
F.write('\tTHETA (deg) \tfi (deg) \tpsi (deg) \tR \tE_MECANICA (J) \tD')

# Las variables que es escriben en el archivo que aparecen a continuación
# se corresponden con el punto final de maniobra del misil para cada
# lanzamiento en cada paso de integración del misil
F.write('\tTime (s)  \tAltura \tDesplazamiento (m)')
F.write('\tVelocidad  \tMach \tFi_l \tPsi_l \tR_l \tEmec Misil \n')

# Cabezas de tabla.
# -------------------------CONDICIONES INICIALES-------------------------
# Ahora, para los próximos cálculos, se definen las variables termodinámicas
# obtenidas del modelo ISA.

# Vectores velocidad, altura, fil y psil vacios, en ellos guardaremos
# los valores de las distintas variables para su posterior uso.

VELOCIDAD = []
ALTURA = []
FI_LIST = []
PSI_LIST = []

Z0 = 12000   # Altitud inicial del misil (m).
R = RT + Z0  # Distancia desde el centro de la Tierra (m).
G0 = gravity(Z0)  # Aceleración gravitatoria inicial (m/s2).
RHO = density(Z0)  # Densidad inicial del aire (kg/m3).
P = pressure(Z0)  # Presión inicial del aire (Pa).
TEMPERATURE = temperature(Z0)  # Temperatura inicial del aire (K).
MU_VISC = viscosity(Z0) # Viscosidad inicial.
W = MASS * G0  # Peso inicial del avión dependiente de la gravedad.

# A la altura inicial el avión vuela en vuelo estacionario.
M1 = 0.9  # Número de Mach inicial.
M = vuelo_crucero(M1)
V = M * (GAMMA * R_AIR * TEMPERATURE)**.5  # Velocidad inicial (m/s).

CL_ALFA1 = cl_alfa(M)  # Pendiente del coeficiente de sustentación.

# Ángulos de asiento, de ataque y de asiento de la velocidad iniciales.

ALFA_NUMERICO = 2 * W / (RHO * V**2 * S_W * CL_ALFA1)
ALFA = angulo_ataque(ALFA_NUMERICO, M)  # Ángulo de ataque.
ALFA_GRADOS = degrees(ALFA)  # Ángulo de ataque en grados.
GAMA = 0  # Ángulo de asiento de la velocidad.
GAMA_GRADOS = degrees(GAMA)  # Ángulo de asiento de la velocidad en grados.
THETA = GAMA + ALFA  # Ángulo de asiento de la velocidad.
THETA_GRADOS = degrees(THETA)  # Ángulo de asiento en grados.

CL = 2 * W / (RHO * V**2 * S_W)
K1 = k(M)
CD01 = cd0(M)
CDMISIL_AVION = cdll(M, V)
CD_INDUCIDA1 = cd_inducida(K1, CL)
CD_INTERFERENCIA1 = CD_interferencia(M)
CD = CD01 + CD_INDUCIDA1  # Polar del avión.  Coeficiente de resistencia.
CD_AVION = CD + 0.35*CD_INTERFERENCIA1    #Añade resistencia de interferencia.


# Fuerzas.

D_AVION = resistencia(V, RHO, CD_AVION)  # Resistencia aerodinámica (N).
D_MISIL = 0.5 * RHO * CDMISIL_AVION * SREF_MISIL * V**2
D = D_AVION + D_MISIL

TH = thrust(M, RHO)  # Empuje (N).
L = sustentacion(V, RHO, CL)

# Esta es la ecuacion en eje horizontal T = D que es la condicion que queremos
# cumplir,por ello calculamos la diferencia y en el while se intenta que sea 0.

DIFERENCIA_TH_D = TH - D

# Ángulos de movimiento bidimensional terrestre iniciales.
PSI = 0      # Ángulo desplazado de la Tierra.
FI = - PSI - GAMA + pi / 2    # Ángulo de la velocidad sobre vertical local.
PSI_GRADOS = degrees(PSI) # Ángulo desplazado de la Tierra en grados.
FI_GRADOS = degrees(FI)   # Ángulo velocidad sobre vertical local en grados.

VELOCIDAD.append(V)
PSI_LIST.append(PSI_GRADOS)
ALTURA.append(Z0)
FI_LIST.append(FI_GRADOS)
# -------------------------INICIO DE LA MANIOBRA-------------------------


# Se impone un valor constante del radio de giro, es decir, giro ascendente
# a radio constante.  En futuras versiones de este cálculo, esto se
# cambiará para buscar una solución más realista.  Por ahora, con objeto de
# facilitar cálculos, se dejará así.
RADIUS = V**2 / (G0 * (N - 1))  # Radio de giro inicial (m).
# Este radio de giro se obtiene para la velocidad inicial en vuelo
# estacionario y para un factor de carga máximo según los pilones de carga
# n = 3,5.
DT = 0.1  # Diferencial de tiempo (s).
'''------SISTEMA DE ECUACIONES PARA PRIMER TRAMO: VUELO ESTACIONARIO------
'''
#Energías.
ECINETICA = .5 * MASS * V**2  # Energía cinética (J).
EPOTENCIAL = MASS * G0 * Z0  # Energía potencial (J).
EMECANICA = ECINETICA + EPOTENCIAL  # Energía mecánica (J).

N = L / W  # Factor de carga.

#Condiciones iniciales para la integración.
T = 0
X = 0
Z = Z0
XT = X
ZT = Z0
OMEGA = V / RADIUS  # Velocidad angular en la maniobra de giro (rad/s).

# Segunda ley de Newton en el eje horizontal (ejes cuerpo).
DV = DT * (TH * cos(ALFA) - D - W * cos(FI)) / MASS

VX = V * sin(FI)  # Velocidad horizontal del avión en ejes HL.
VY = V * cos(FI)  # Velocidad vertical del avión en ejes VL.

VXT = V * cos(GAMA) # Velocidad horizontal del avión en ejes tierra (m/s).
VYT = V * sin(GAMA) # Velocidad vertical del avión en ejes tierra (m/s).

DX = VX * DT  # Variación horizontal de la posición en "Horizonte local".
DZ = VY * DT  # Variación vertical de la posición en "Vertical local".

DXT = VXT * DT  # Variación horizontal de la posición en ejes tierra (m/s).
DZT = VYT * DT  # Variación vertical de la posición en ejes tierra (m/s).

DTHETA = V * DT / RADIUS  # Variación del ángulo de asiento.
DPSI = ((V*sin(FI))/(R))*DT # Variación del ángulo de la vertical local.

# Incialización de las primeras variables del misil.
VL = V    # Se declara que la velocidad de lanz. del misil es la del avión.
ML = VL / ((GAMMA * R_AIR * TEMPERATURE)**0.5) # Mach para velocidad del misil.
CDL = cdll(ML, Z)  # Cd del misil para la cond. de vuelo.
D_MISIL = 0.5 * RHO * CDL * SREF_MISIL * VL**2  # Resis. del misil.


# -------SISTEMA DE ECUACIONES PARA SEGUNDO TRAMO: MANIOBRA DE GIRO-------
# Ahora comienza el bucle relativo al giro ascendente, que analiza la
# trayectoria con nuevas ecuaciones y condiciones de vuelo que se detallan
# más afelante.  El significado de theta < beta implica que el bucle realice
# el cálculo requerido siempre que el ángulo theta sea menor que beta.  Se ha
# obligado a que beta sea el ángulo de final de giro (al inicio del programa
# se le ha dado un rango de valores).  Por tanto, una vez que theta sea igual
# a beta, se dará por concluida la maniobra de giro y se comenzará con el
# siguiente tramo. Todo ello mientras la velocidad sea positiva, condición
# que permite ver si el resultado es realista o no.

while GAMA < BETA and V > 0:
    # A continuación, se escriben en un fichero todas las variables aquí
    # detalladas para cada valor de theta < beta.

    F.write('%.8f\t' %T)  # Tiempo (s).
    F.write('%.8f\t' %Z)  # Altitud (m).
    F.write('%.8f\t' %X) # Recorrido (m).
    F.write('%.8f\t' %V) # Velocidad (m/s).
    F.write('%.8f\t' %M)  # Número de Mach.
    F.write('%.8f\t' %CD_AVION)    # Coeficiente de Resistencia.
    F.write('%.8f\t' %D)
    F.write('%.8f\t' %TH)
    F.write('%.8f\t' %ALFA_GRADOS)  # Ángulo de ataque (grados).
    F.write('%.8f\t' %GAMA_GRADOS)  # Asiento de la velocidad (grados).
    F.write('%.8f\t' %THETA_GRADOS)  # Ángulo de asiento (grados).
    F.write('%.8f\t' %FI_GRADOS)     # Ángulo del avión sobre la VL.
    F.write('%.8f\n' %PSI_GRADOS)    # Ángulo de la vertical local.

    # Ya que este análisis de maniobra, a diferencia del anterior, lleva un
    # cálculo para distintos valores de tiempo y velocidad, se debe
    # programar su evolución en términos de sus variaciones diferenciales.
    # Aquí se escriben los incrementos diferenciales de las coordenadas
    # espaciales, temporales y de velocidad.  Los diferenciales se obtienen
    # del tramo anterior, y sirven para que estos nuevos valores actúen en
    # las nuevas condiciones para calcular nuevas variables.

    # Inicio del cálculo de la maniobra del avión.

    # Incrementos.

    T = T + DT  # Evolución temporal (s).
    V = V + DV  # Velocidad (m/s).
    X = X + DX  # Posición horizontal-local (m).
    XT = XT + DXT  # Posición horizontal en ejes tierra (m).
    Z = Z + DZ  # Altitud en horizonte-local (m).
    ZT = ZT + DZT  # Altitud en ejes del punto de lanzamiento.

    # Guardamos las variables en la lista.
    VELOCIDAD.append(V)
    ALTURA.append(Z)

    # Gravitación
    R = RT + Z  # Distancia al centro de la Tierra (m).
    G = gravity(Z)  # Aceleración de la gravedad (m/s2).
    W = MASS * G    # Peso dependiente de la gravedad

    # Las variables termodinámicas habrán variado con la altura.
    RHO = density(Z)  # Densidad (kg/m3).
    TEMPERATURE = temperature(Z)  # Temperatura (K).
    MU_VISC = viscosity(Z) # Viscosidad

    M = V / (GAMMA * R_AIR * TEMPERATURE)**.5 # Mach de vuelo.
    N = 3.5  # Tomamos la condición de factor de carga máximo y constante.
    RADIUS = V**2 / (G * (N - 1))  # Radio de giro varía con la velocidad.

    #Las características aerodinámicas varían con el número de Mach.
    CL_ALFA1 = cl_alfa(M)
    K1 = k(M)
    CD01 = cd0(M)
    ALFA_NUMERICO = 2 * W * N / (RHO * V**2 * S_W * CL_ALFA1)
    #El nuevo ángulo de ataque resultará del nuevo factor de carga (antes 1
    # y ahora máximo), la nueva velocidad y las nuevas características
    # aerodinámicas.
    ALFA = angulo_ataque(ALFA_NUMERICO, M)  # Ángulo de ataque.
    THETA = THETA + DTHETA  # Ángulo de asiento (empuje horizontal).
    GAMA = THETA - ALFA  # Ángulo de asiento de la velocidad.
    PSI = PSI + DPSI     # Ángulo de vertical local
    FI = - PSI - GAMA + pi / 2  # Ángulo de la vel sobre la vertical local

    # Ángulos en grados para la exportación a los ficheros.
    ALFA_GRADOS = degrees(ALFA)  # Ángulo de ataque (grados).
    THETA_GRADOS = degrees(THETA)  # Ángulo de asiento (grados).
    GAMA_GRADOS = degrees(GAMA)  # Asiento de la velocidad (grados).
    PSI_GRADOS = degrees(PSI)    # Ángulo de vertical local (grados).
    FI_GRADOS = degrees(FI)  # Ángulo velocidad sobre la vert loca(grados).

    'Guardamos las variables en la lista'
    PSI_LIST.append(PSI_GRADOS)
    FI_LIST.append(FI_GRADOS)

    CL = ALFA * CL_ALFA1  # Coeficiente de sustentación.
    N = .5 * RHO * V**2 * S_W * CL / W

    # Nuevas proyecciones de la velocidad
    # (ejes radiales y ejes posición inicial).

    VX = V * sin(FI) # Velocidad horizontal en ejes locales
    VY = V * cos(FI) # Velocidad vertical en ejes locales
    VXT = V * cos(GAMA) # Proyección horizontal velocidad en ejes tierra
    VYT = V * sin(GAMA) # Proyección vertical velocidad en ejes tierra
    CD_INDUCIDA1 = cd_inducida(K1, CL) # Coef. resistencia inducida
    CD_INTERFERENCIA1 = CD_interferencia(M)
    CD = CD01 + CD_INDUCIDA1  # Coeficiente de resistencia total
    CD_AVION = CD + 0.35 * CD_INTERFERENCIA1
    D = resistencia(V, RHO, CD_AVION)  # Fuerza de resistencia (N).
    L = N * W  # Fuerza de sustentación (N).

    #Energías.
    ECINETICA = .5 * MASS * V**2  # Energía cinética (J).
    EPOTENCIAL = MASS * G * Z   # Energía potencial (J).
    EMECANICA = ECINETICA + EPOTENCIAL  # Energía mecánica (J).

    #Empuje.
    TH = thrust(M, RHO)  # Empuje (N).
    DIFERENCIA_TH_D = TH - D # Diferencia de empuje frente a resitencia.

    #Velocidad angular.
    OMEGA = V / RADIUS  # Velocidad angular (rad/s).

    #Nuevas variaciones diferenciales.  Segunda ley de Newton.
    DV = DT * (TH * cos(ALFA) - D - W * cos(FI)) / MASS
    DPSI = ((V * sin(FI)) / (R)) * DT
    #Esta última ecuación nos permite obtener un nuevo incvremento de
    # velocidad (m/s) que, al reiniciar el bucle será sumado al último
    # valor de la velocidad, obteniendo la variación de la velocidad.

    # Diferenciales
    DX = VX * DT  # Variación horizontal de la posición (m).
    DZ = VY * DT  # Variación vertical de la posición (m).

    DXT = VXT * DT # Variación horizontal posición en ejes tierra (m/s).
    DZT = VYT * DT # Variación vertical posición en ejes tierra (m/s).

    DTHETA = OMEGA * DT  # Variación del ángulo de asiento.


F.close()

# Como resumen:
# 1) El código ha empezado en una condición de vuelo uniforme.
# 2) La siguiente maniobra es un giro ascendente, a factor de carga máximo y
# constante, y con mínima resistencia.
# 3) La última maniobra es un ascenso con el ángulo final del giro, con
# coeficiente de sustentación óptimo.
# Este programa nos exportará un total de 9 archivos (uno para cada incremento
# de 10º del ángulo final de giro) que, exportados a Excel nos permiten
# observar cómo cambian las variables según las condiciones de vuelo.
