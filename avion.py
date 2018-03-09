# -*- coding: utf-8 -*-
"""
@author: Team REOS
"""
from math import radians, cos, sin, degrees, pi
from modeloISA import density, temperature, GAMMA, viscosity, pressure, R_AIR, gravity
from modelo_empuje import thrust
from aero_avion import cl_alfa, angulo_ataque, k, cd0, cd_inducida, S_W
from aero_avion import resistencia, sustentacion
from aero_misil import cdll, SREF_MISIL, SGASES
from velocidad_crucero import vuelo_crucero


#--------------------------CONDICIONES GRAVITATORIAS--------------------------
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2).


#------------------------CARACTERÍSTICAS DE LA AERONAVE------------------------

N = 3.5  # Factor de carga máximo.
MASS = 14273  # Masa del avión cargado (kg).
W = MASS * GRAV  # Peso del avión (N).

beta=89 
                       
beta = radians(beta)
f = open('resultados', 'w')  # Fichero de escritura sin extensión.

'''Las variables que se escriben a continuación corresponden al paso de 
integración del avión en cada paso'''
f.write('TIME (s)\tALTURA lanzamiento (m) \tDesplazamiento (m)')
f.write ('\tVELOCIDAD (m/s)\tMACH\tALFA (deg) \tGAMMA (deg) ')
f.write('\tTHETA (deg) \tfi (deg) \tpsi (deg) \tR \tE_MECANICA (J) \tD')

'''Las variables que es escriben en el archivo que aparecen a continuación 
se corresponden con el punto final de maniobra del misil para cada 
lanzamiento en cada paso de integración del misil'''
f.write('\tTime (s)  \tAltura \tDesplazamiento (m)')
f.write('\tVelocidad  \tMach \tFi_l \tPsi_l \tR_l \tEmec Misil \n')

# Cabezas de tabla.
'''-------------------------CONDICIONES INICIALES-------------------------
Ahora, para los próximos cálculos, se definen las variables termodinámicas
obtenidas del modelo ISA.
'''
'''Vectores velocidad, altura, fil y psil vacios, en ellos guardaremos los valores
de las distintas variables para su posterior uso'''

velocidad = []
altura = []
fi_list = []
psi_list = []

z0 = 12000   # Altitud inicial del misil (m).
r = RT + z0  # Distancia desde el centro de la Tierra (m).
g0 = gravity(z0)  # Aceleración gravitatoria inicial (m/s2).
rho = density(z0)  # Densidad inicial del aire (kg/m3).
p = pressure(z0)  # Presión inicial del aire (Pa).
T = temperature(z0)  # Temperatura inicial del aire (K).
Mu_Visc = viscosity(z0) # Viscosidad inicial 
W = MASS * g0  # Peso inicial del avión dependiente de la gravedad

#A la altura inicial el avión vuela en vuelo estacionario.
M1 = 1.4  # Número de Mach inicial.
M = vuelo_crucero(M1)
v = M * (GAMMA * R_AIR * T)**.5  # Velocidad inicial (m/s).

CL_alfa1 = cl_alfa(M)  # Pendiente del coeficiente de sustentación.

# Ángulos de asiento, de ataque y de asiento de la velocidad iniciales.
alfa_numerico = 2 * W / (rho * v**2 * S_W * CL_alfa1)
alfa = angulo_ataque(alfa_numerico, M)  # Ángulo de ataque.
alfa_grados = degrees(alfa)  # Ángulo de ataque en grados.
gama = 0  # Ángulo de asiento de la velocidad.
gama_grados = degrees(gama)  # Ángulo de asiento de la velocidad en grados.
theta = gama + alfa  # Ángulo de asiento de la velocidad.
theta_grados = degrees(theta)  # Ángulo de asiento en grados.

CL = 2 * W / (rho * v**2 * S_W)
k1 = k(M)
CD01 = cd0(M)
CDmisilavion = cdll(M, v)
CD_inducida1 = cd_inducida(k1, CL)
CD = CD01 + CD_inducida1  # Polar del avión.  Coeficiente de resistencia.


#Fuerzas.
Davion = resistencia(v, rho, CD)  # Resistencia aerodinámica (N).
Dmisil = 0.5 * rho * CDmisilavion * SREF_MISIL * v**2
D = Davion + Dmisil

Th = thrust(M, rho)  # Empuje (N).

''' Esta es la ecuacion en eje horizontal T = D que es la condicion que queremos cumplir
por ello calculamos la diferencia y en el while se intenta que sea 0 '''

diferencia_T_D = Th - D

# Ángulos de movimiento bidimensional terrestre iniciales
psi = 0      # Ángulo desplazado de la Tierra
fi = - psi - gama + pi / 2    # Ángulo de la velocidad sobre vertical local
psi_grados = degrees(psi) # Ángulo desplazado de la Tierra en grados
fi_grados = degrees(fi)   # Ángulo velocidad sobre vertical local en grados

velocidad.append(v)
psi_list.append(psi_grados)
altura.append(z0)
fi_list.append(fi_grados)
'''-------------------------INICIO DE LA MANIOBRA-------------------------
'''

#Se impone un valor constante del radio de giro, es decir, giro ascendente
# a radio constante.  En futuras versiones de este cálculo, esto se
# cambiará para buscar una solución más realista.  Por ahora, con objeto de
# facilitar cálculos, se dejará así.
radius = v**2 / (g0 * (N - 1))  #Radio de giro inicial (m).
# Este radio de giro se obtiene para la velocidad inicial en vuelo
# estacionario y para un factor de carga máximo según los pilones de carga
# n = 3,5.
dt = 0.1  # Diferencial de tiempo (s).
'''------SISTEMA DE ECUACIONES PARA PRIMER TRAMO: VUELO ESTACIONARIO------
'''
#Energías.
ecinetica = .5 * MASS * v**2  # Energía cinética (J).
epotencial = MASS * g0 * z0  # Energía potencial (J).
emecanica = ecinetica + epotencial  # Energía mecánica (J).

n = L / W  # Factor de carga.

#Condiciones iniciales para la integración.
t = 0
x = 0
z = z0
xt = x
zt = z0
omega = v / radius  # Velocidad angular en la maniobra de giro (rad/s).

#Segunda ley de Newton en el eje horizontal (ejes cuerpo).
dv = dt * (Th * cos(alfa) - D - W * cos(fi)) / MASS

vx = v * sin(fi)  # Velocidad horizontal del avión en ejes HL.
vy = v * cos(fi)  # Velocidad vertical del avión en ejes VL.

vxt = v * cos(gama) # Velocidad horizontal del avión en ejes tierra (m/s).
vyt = v * sin(gama) # Velocidad vertical del avión en ejes tierra (m/s).

dx = vx * dt  # Variación horizontal de la posición en "Horizonte local".
dz = vy * dt  # Variación vertical de la posición en "Vertical local".

dxt = vxt * dt  # Variación horizontal de la posición en ejes tierra (m/s).
dzt = vyt * dt  # Variación vertical de la posición en ejes tierra (m/s).

dtheta = v * dt / radius  # Variación del ángulo de asiento.
dpsi = ((v*sin(fi))/(r))*dt # Variación del ángulo de la vertical local.

# Incialización de las primeras variables del misil   
vl = v    # Se declara que la velocidad de lanz. del misil es la del avión.
Ml = vl / ((GAMMA * R_AIR * T)**0.5) # Mach para dicha velocidad del misil.
Cdl = cdll(Ml, z)  # Cd del misil para la cond. de vuelo.
D_misil = 0.5 * rho * Cdl * SREF_MISIL * vl**2  # Resis. del misil.


'''-------SISTEMA DE ECUACIONES PARA SEGUNDO TRAMO: MANIOBRA DE GIRO-------
Ahora comienza el bucle relativo al giro ascendente, que analiza la
trayectoria con nuevas ecuaciones y condiciones de vuelo que se detallan
más afelante.  El significado de theta < beta implica que el bucle realice
el cálculo requerido siempre que el ángulo theta sea menor que beta.  Se ha
obligado a que beta sea el ángulo de final de giro (al inicio del programa
se le ha dado un rango de valores).  Por tanto, una vez que theta sea igual
a beta, se dará por concluida la maniobra de giro y se comenzará con el
siguiente tramo. Todo ello mientras la velocidad sea positiva, condición
que permite ver si el resultado es realista o no.
'''

while gama < beta and v > 0:
    '''A continuación, se escriben en un fichero todas las variables aquí
    detalladas para cada valor de theta < beta.
    '''
    f.write('%.8f\t' %t)  # Tiempo (s).
    f.write('%.8f\t' %z)  # Altitud (m).
    f.write('%.8f\t' %x)  # Recorrido (m)
    f.write('%.8f\t' %v)  # Velocidad (m/s).
    f.write('%.8f\t' %M)  # Número de Mach.
    f.write('%.8f\t' %alfa_grados)  # Ángulo de ataque (grados).
    f.write('%.8f\t' %gama_grados)  # Asiento de la velocidad (grados).
    f.write('%.8f\t' %theta_grados)  # Ángulo de asiento (grados).
    f.write('%.8f\t' %fi_grados)     # Ángulo del avión sobre la VL.
    f.write('%.8f\n' %psi_grados)    # Ángulo de la vertical local.
    #Ya que este análisis de maniobra, a diferencia del anterior, lleva un
    # cálculo para distintos valores de tiempo y velocidad, se debe
    # programar su evolución en términos de sus variaciones diferenciales.
    #Aquí se escriben los incrementos diferenciales de las coordenadas
    # espaciales, temporales y de velocidad.  Los diferenciales se obtienen
    # del tramo anterior, y sirven para que estos nuevos valores actúen en
    # las nuevas condiciones para calcular nuevas variables.
   
    
    '''Inicio del cálculo de la maniobra del avión '''
    
    #Incrementos. 
    t = t + dt  # Evolución temporal (s).
    v = v + dv  # Velocidad (m/s).
    x = x + dx  # Posición horizontal-local (m).
    xt = xt + dxt  # Posición horizontal en ejes tierra (m)
    z = z + dz  # Altitud en horizonte-local (m).
    zt = zt + dzt  # Altitud en ejes del punto de lanzamiento.
    
    'Guardamos las variables en la lista'
    velocidad.append(v)
    altura.append(z)
    
    # Gravitación
    r = RT + z  # Distancia al centro de la Tierra (m).
    g = gravity(z)  # Aceleración de la gravedad (m/s2).
    W = MASS * g    # Peso dependiente de la gravedad
    
    # Las variables termodinámicas habrán variado con la altura.
    rho = density(z)  # Densidad (kg/m3).
    T = temperature(z)  # Temperatura (K).
    Mu_Visc = viscosity(z) # Viscosidad
    
    M = v / (GAMMA * R_AIR * T)**.5 # Mach de vuelo.
    n = 3.5  # Tomamos la condición de factor de carga máximo y constante.
    radius = v**2 / (g * (n - 1))  # Radio de giro varía con la velocidad.
    
    #Las características aerodinámicas varían con el número de Mach.
    CL_alfa1 = cl_alfa(M)
    k1 = k(M)
    CD01 = cd0(M)
    alfa_numerico = 2 * W * n / (rho * v**2 * S_W * CL_alfa1)
    #El nuevo ángulo de ataque resultará del nuevo factor de carga (antes 1
    # y ahora máximo), la nueva velocidad y las nuevas características
    # aerodinámicas.
    alfa = angulo_ataque(alfa_numerico, M)  # Ángulo de ataque.
    theta = theta + dtheta  # Ángulo de asiento (empuje horizontal).
    gama = theta - alfa  # Ángulo de asiento de la velocidad.
    psi = psi + dpsi     # Ángulo de vertical local
    fi = - psi - gama + pi / 2  # Ángulo de la vel sobre la vertical local
    
    # Ángulos en grados para la exportación a los ficheros.
    alfa_grados = degrees(alfa)  # Ángulo de ataque (grados).
    theta_grados = degrees(theta)  # Ángulo de asiento (grados).
    gama_grados = degrees(gama)  # Asiento de la velocidad (grados).
    psi_grados = degrees(psi)    # Ángulo de vertical local (grados).
    fi_grados = degrees(fi)  # Ángulo velocidad sobre la vert loca(grados).
    
    'Guardamos las variables en la lista'
    psi_list.append(psi_grados)
    fi_list.append(fi_grados)
    
    CL = alfa * CL_alfa1  # Coeficiente de sustentación.
    n = .5 * rho * v**2 * S_W * CL / W
    
    #Nuevas proyecciones de la velocidad (ejes radiales y ejes posición inicial).
    
    vx = v * sin(fi) # Velocidad horizontal en ejes locales                                                            
    vy = v * cos(fi) # Velocidad vertical en ejes locales
    vxt = v * cos(gama) # Proyección horizontal velocidad en ejes tierra
    vyt = v * sin(gama) # Proyección vertical velocidad en ejes tierra
    CD_inducida1 = cd_inducida(k1, CL) # Coef. resistencia inducida
    CD = CD01 + CD_inducida1  # Coeficiente de resistencia total
    D = resistencia(v, rho, CD)  # Fuerza de resistencia (N).
    L = n * W  # Fuerza de sustentación (N).
    
    #Energías.
    ecinetica = .5 * MASS * v**2  # Energía cinética (J).
    epotencial = MASS * g * z   # Energía potencial (J).
    emecanica = ecinetica + epotencial  # Energía mecánica (J).
    
    #Empuje.
    Th = thrust(M, rho)  # Empuje (N).
    diferencia_T_D = Th - D # Diferencia de empuje frente a resitencia. 
    
    #Velocidad angular.
    omega = v / radius  # Velocidad angular (rad/s).
    
    #Nuevas variaciones diferenciales.  Segunda ley de Newton.
    dv = dt * (Th * cos(alfa) - D - W * cos(fi)) / MASS
    dpsi = ((v * sin(fi)) / (r)) * dt
    #Esta última ecuación nos permite obtener un nuevo incvremento de
    # velocidad (m/s) que, al reiniciar el bucle será sumado al último
    # valor de la velocidad, obteniendo la variación de la velocidad.
    
    # Diferenciales
    dx = vx * dt  # Variación horizontal de la posición (m).
    dz = vy * dt  # Variación vertical de la posición (m).
    
    dxt = vxt * dt # Variación horizontal posición en ejes tierra (m/s).
    dzt = vyt * dt # Variación vertical posición en ejes tierra (m/s).
    
    dtheta = omega * dt  # Variación del ángulo de asiento       
    

f.close()
           
'''
Como resumen:
1) El código ha empezado en una condición de vuelo uniforme.
2) La siguiente maniobra es un giro ascendente, a factor de carga máximo y
constante, y con mínima resistencia.
3) La última maniobra es un ascenso con el ángulo final del giro, con
coeficiente de sustentación óptimo.

Este programa nos exportará un total de 9 archivos (uno para cada incremento de
10º del ángulo final de giro) que, exportados a Excel nos permiten observar
cómo cambian las variables según las condiciones de vuelo.
'''
    

            
        
