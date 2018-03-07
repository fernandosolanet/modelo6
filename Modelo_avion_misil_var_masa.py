# -*- coding: utf-8 -*-
"""
@author: Team REOS

Este programa servirá para establecer la trayectoria a describir por el
avión y la puesta en órbita del cohete REOS.  Contiene el loop del avión y el 
lanzamiento del cohete de una sola etapa.  Se extraen en un archivo de texto
las variables de interés durante el vuelo para su posterior análisis.

"""

from math import radians, cos, pi, sin, degrees, atan, log10
from modeloISA import density, temperature, GAMMA, viscosity, pressure, R_AIR
from modelo_empuje import thrust
from aero_avion import cl_alfa, angulo_ataque, k, cd0, cd_inducida
from aero_avion import resistencia, sustentacion


#--------------------------CONDICIONES GRAVITATORIAS--------------------------
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2).


#-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL VEHÍCULO-------------------


N = 3.5  # Factor de carga máximo.
S_W = 49.2  # Superficie alar (m2).
B = 11.7  # Envergadura (m).
AR = B**2 / S_W  # Alargamiento = 2,78.
FLECHA = radians(52)  # Flecha BA.
FLECHA2 = radians(41.4)  # Flecha 1/4.
LF = 19.2  # Longitud del fuselaje (m).
BF = 2.87  # Longitud del fuselaje (m).
#Se desprecian las pérdidas por consumo de combustible en el peso del avión.

K = .4  # EL perfil del F-4 es el NACA0006.4-64 por tanto la K es 0,4.
ESTRECHAMIENTO = .26  # Estrechamiento.
ESPESOR = .064  # Espesor relativo máximo.

NE = 2  # Número de motores.
DM = 0  # Diámetro de cada motor.
LM = 0  # Longitud de cada motor.
#La longitud y el diámetro de cada motor aparecen nulos porque no se ha
# refelejado todavía en los cálculos la importancia de la geometría de los
# motores.

MASS = 14273  # Masa del avión cargado (kg).
W = MASS * GRAV  # Peso del avión (N).

#A continuación, se definen las superficies de mando del avión, que nos
# servirán para, más adelante, calcular el coeficiente de resistencia parásita.
S_LEX = 0  # Área del Lex. [¡!]
S_H = 6.39  # Área de la superficie de mando horizontal (m2).
S_V = 5.035  # Área de la superficie de mando vertical (m2).

#Se necesitarán para más adelante los valores de Mach crítico y de divergencia,
# los cuales son función de la flecha y del espesor relativo máximo.  Estos
# valores marcarán los límites de los dominios que nos servirán para calcular
# el coeficicente de resitencia inducida con efecto de compresibilidad.
M_C = 1 - (.065 * (cos(FLECHA))**.6) * (100 * ESPESOR)**.637814
# Mach crítico.
M_D = 1.08 * M_C  # Mach de divergencia.
M_D098 = .98 * M_D


'''
-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL MISIL-------------------
'''

DIAMETRO_M = .5  # Diámetro del misil (m).
LONGITUD_CONO = .9  # Longitud del cono del misil (m).
LONGITUD_MISIL = 3  # Longitud total del misil (m).

SW_ALETA = .07875  # Superficie de una aleta del AIM (tomado como ref., m2).
CRAIZ_ALETA = .24  # Cuerda raiz de la aleta (m).
CMEDIA_ALETA = .18  # Cuerda media de la aleta (m).
ESPESOR_ALETA = .0065  # Espesor de la aleta (m).
TAO_ALETA = ESPESOR_ALETA / CMEDIA_ALETA  # TAO de la aleta.
NUM_ALETAS = 4  # Número de aletas.
SWTOTAL_ALETAS = SW_ALETA * NUM_ALETAS  # Superficie total de aletas (m2).

SUP_CONO = pi * DIAMETRO_M / 2 * (LONGITUD_CONO**2 + DIAMETRO_M**2 / 4)**.5
# Superficie exterior del cono (m2).
SUP_TOTAL = pi * DIAMETRO_M * (LONGITUD_MISIL - LONGITUD_CONO)
# Superficie exterior del misil (m2).
SREF_MISIL = pi * DIAMETRO_M**2 / 4  # Superficie de referencia del misil (m2).
SGASES = pi * (DIAMETRO_M * .45)**2
# Área de salida de los gases (consideramos el área de salida de la tobera,
# m2).
RATIO_AREAS = 1 - SGASES / SREF_MISIL  # Relación de áreas.
ANGULO_CONO = degrees(atan(.5 * DIAMETRO_M / LONGITUD_CONO))
# Ángulo del cono (deg).


#----------------------------AERODINÁMICA DEL MISIL----------------------------

def coef_resistencia_base_misil(mach):
    '''Coeficiente de resistencia base del misil.  Varía con el número de Mach.
    '''
    if mach < .8:
        return 0
    elif mach < 1:
        x_0 = -1.548523
        x_1 = 6.05972764
        x_2 = -7.30548391
        x_3 = 2.96129532
        x_4 = 0
    elif mach < 1.1:
        x_0 = 5790.90984
        x_1 = -21984.3314
        x_2 = 31277.4812
        x_3 = -19764.4892
        x_4 = 4680.59822
    elif mach < 1.5:
        x_0 = -4.11856506
        x_1 = 14.2267421
        x_2 = -16.9678524
        x_3 = 8.771665
        x_4 = -1.67398037
    elif mach < 2.2:
        x_0 = .30748
        x_1 = -.13258
        x_2 = .028812
        x_3 = 0
        x_4 = 0
    elif mach > 2.2:
        x_0 = .18481
        x_1 = -.022895
        x_2 = .0051876
        x_3 = -.00040742
        x_4 = 0
    elif mach > 3.5:
        return .15
    return x_4 * mach**4 + x_3 * mach**3 + x_2 * mach**2 + x_1 * mach + x_0

def cfcono_misil(re_cono, machl):
    '''Coeficiente de fricción del cono.
    '''
    #LAMINAR
    if re_cono < 1e6:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfi_cono = .664 * re_cono**(-1 / 2)
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO.
        cf_cono = 2 * cfi_cono
        #CÁLCULO COEFICIENTE DE FRICCIÓN COMPRESIBLE.
        cfm_cono = cf_cono / (1 + .17 * machl**2)**.1295
    #TURBULENTO
    else:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfi_cono = .288 / log10(re_cono)**2.45
        #CALCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE.
        cf_cono = cfi_cono * 1.597 / log10(re_cono)**.15
        #CÁLCULO COEFICIENTE DE FRICCIÓN MEDIO.
        cfm_cono = cf_cono / (1 + (GAMMA - 1) / 2 * machl**2)**.467
    return cfm_cono * SUP_CONO / SREF_MISIL

def cfcil(re_cilindro, machl):
    '''Coeficiente de fricción del cilindro.
    '''
    #LAMINAR
    if re_cilindro < 1e6:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfi_cil = .664 * re_cilindro**(-1 / 2)
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO.
        cf_cil = 2 * cfi_cil
        #CÁLCULO COEFICIENTE DE FRICCIÓN COMPRESIBLE.
        cfm_cil = cf_cil / (1 + .17 * machl**2)**.1295
    #TURBULENTO
    else:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfi_cil = .288 / log10(re_cilindro)**2.45
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE.
        cf_cil = cfi_cil * 1.597 / log10(re_cilindro)**.15
        #CÁLCULO COEFICIENTE DE FRICCIÓN MEDIO.
        cfm_cil = cf_cil / (1 + (GAMMA - 1) / 2 * machl**2)**.467
    return cfm_cil * SUP_TOTAL / SREF_MISIL

def cd_wave(mach, angulo, cd_f):
    '''Coeficiente de onda del misil.
    '''
    if mach >= 1:
        #RÉGIMEN SUPERSÓNICO.
        return (.083 + .096 / mach**2) * (angulo / 10)**1.69
    #RÉGIMEN SUBSÓNICO.
    ratio = LONGITUD_CONO / DIAMETRO_M
    return (60 / ratio**3 + .0025 * ratio) * cd_f

def cd_wave_aletas(mach):
    '''Coeficiente de onda de las aletas.
    '''
    #RÉGIMEN SUPERSÓNICO
    if mach >= 1:
        return 4 * TAO_ALETA**2 / (mach**2 - 1)**.5 * (SWTOTAL_ALETAS
                                                       / SREF_MISIL)
    #RÉGIMEN SUBSÓNICO.
    return 0

def cf_aletas(reyn_aleta, mach):
    '''Coeficiente de fricción de las aletas.
    '''
    #LAMINAR.
    if reyn_aleta < 1e6:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfialetas = .664 / reyn_aleta**.5
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO.
        cf1aletas = 2 * cfialetas
        #CÁLCULO COEFICIENTE DE FRICCIÓN COMPRESIBLE.
        cfmaletas = cf1aletas / (1 + .17 * mach**2)**.1295
    #TURBULENTO.
    else:
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfialetas = .288 * (log10(reyn_aleta))**(-2.45)
        #CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE.
        cf1aletas = cfialetas * 1.597 * ((log10(reyn_aleta))**(-.15))
        #CÁLCULO COEFICIENTE DE FRICCIÓN MEDIO.
        cfmaletas = cf1aletas / (1 + (GAMMA - 1) / 2 * mach**2)**.467
    return cfmaletas * SWTOTAL_ALETAS / SREF_MISIL

def cdll(machl):
    '''Coeficiente de resistencia total del misil.
    '''
    #CÁLCULO DEL COEFICIENTE DE RESISTENCIA BASE.
    cd_base_misil = coef_resistencia_base_misil(machl) * RATIO_AREAS
    re_cono = rho * vl * LONGITUD_CONO / Mu_Visc
    # Número de Reynolds en el cono.
    re_cil = rho * vl * (LONGITUD_MISIL - LONGITUD_CONO) / Mu_Visc
    # Número de Reynolds en el cilindro.
    #CÁLCULO DEL COEFICIENTE DE FRICCIÓN TOTAL REFERIDO A LA SUPERFICIE
    # TRANSVERSAL.
    cd_friccion_cono = cfcono_misil(re_cono, machl)
    # Coeficiente de fricción en el cono.
    cd_friccion_cil = cfcil(re_cil, machl)
    # Coeficiente de fricción en el cilindro.
    cd_friccion = cd_friccion_cono + cd_friccion_cil
    #CÁLCULO DEL COEFICIENTE DE ONDA.
    cd_onda = cd_wave(machl, ANGULO_CONO, cd_friccion)
    #RESISTENCIA DE LAS ALETAS.
    ##COEFICIENTE DE ONDA.
    cd_onda_aletas = cd_wave_aletas(machl)
    re_aletas = rho * vl * CRAIZ_ALETA / Mu_Visc
    # Número de Reynolds en las aletas.
    #COEFICIENTE DE FRICCIÓN DE LAS ALETAS.
    cdfriccion_aletas = cf_aletas(re_aletas, machl)
    return (cd_base_misil + cd_friccion + cd_onda + cd_onda_aletas
            + cdfriccion_aletas)

#-------------------------TRAYECTORIA/PUESTA EN ÓRBITA-------------------------
#En esta parte del código se escribe la integración de las ecuaciones de los
# distintos tramos del vuelo y la puesta en órbita del cohete.  Se comienza
# abriendo un bucle para distintos valores del impulso específico del
# propulsante del cohete REOS.  Este bucle debe durar todo el programa para
# conseguir que exporte los distintos ficheros correspondientes a cada Isp.

for isp in range(int(220 * 9.8), int(320 * 9.8), 100):
    #Bucle para el impulso específico, desde un valor en segundos de 220
    # hasta un valor de 320. Así se verá la influencia de este parámetro
    # en la puesta en órbita del cohete.
    gasto = 60  # Gasto másico del misil (kg/s).
    #Masa total del misil: 1000 kg.  Definido más adelante.
    masa_propulsante = 750  # Masa de propulsante (kg).
    empuje_misil = gasto * isp
    # Empuje variable para cada ensayo (varía con el Isp a gasto cte).
    t_combustion = masa_propulsante / gasto  # Tiempo de combustión (s).
    beta = 89
    #El ángulo beta determina el ángulo final de la maniobra de giro; es
    # decir, es el ángulo de asiento del avión con el que se quiere que
    # comience el ascenso tras el giro en el plano vertical.
    beta = radians(beta)
    isp_texto = str(isp)
    # Para el nombre del archivo de texto que recoge los datos
    f = open(isp_texto, 'w')  # Fichero de escritura sin extensión.
    f.write('TIEMPO DE LANZAMIENTO(s)\tALTURA DE LANZAMIENTO (m)\tVELOCIDAD')
    f.write('DE LANZAMIENTO(m/s)\tMACH\tALFA (deg)\tGAMMA (deg)\tTHETA (deg)')
    f.write('\tE_MECÁNICA (J)\tRESISTENCIA (N)\tTIEMPO (s)\tALTURA (m)\t')
    f.write('VELOCIDAD (m/s)\tMACH\tTHETA (deg)\tMASA (kg)\tE_MECÁNICA (J)\t')
    f.write('POSICIÓN HORIZONTAL FINAL (m)\tEMPUJE (N)\tIMPULSO ESPECÍFICO ')
    f.write('(N s/kg)\n')  # Cabezas de tabla.
    #--------------------------CONDICIONES INICIALES---------------------------
    #Ahora, para los próximos cálculos, se definen las condiciones iniciales de
    # las variables.
    h = 12000  # Altitud inicial (m).
    r = RT + h  # Distancia desde el centro de la Tierra (m).
    g0 = MU / r**2  # Aceleración gravitatoria (m/s2).
    rho = density(h)  # Densidad inicial del aire (kg/m3).
    p = pressure(h)  # Presión inicial del aire (Pa).
    T = temperature(h)  # Temperatura inicial del aire (K).
    Mu_Visc = viscosity(h)  # Viscosidad (Pa s).
    #A la altura inicial el avión vuela en vuelo estacionario.
    M = 1.8  # Número de Mach inicial.
    v = M * (GAMMA * R_AIR * T)**.5  # Velocidad inicial (m/s).
    CL_alfa1 = cl_alfa(M)  # Pendiente del coeficiente de sustentación.
    #Ángulos de asiento, de ataque y de asiento de la velocidad iniciales.
    alfa_numerico = 2 * W / (rho * v**2 * S_W * CL_alfa1)
    alfa = angulo_ataque(alfa_numerico, M)  # Ángulo de ataque.
    alfa_grados = degrees(alfa)  # Ángulo de ataque (deg).
    gama = 0  # Ángulo de asiento.
    gama_grados = degrees(gama)  # Ángulo de asiento (deg).
    theta = gama + alfa  # Ángulo de asiento de la velocidad.
    theta_grados = degrees(theta)  # Ángulo de asiento de la velocidad (deg).
    #Coeficientes aerodinámicos.
    CL = alfa * CL_alfa1  # Coeficiente de sustentación inicial.
    k1 = k(M)
    CD01 = cd0(M)
    CD_inducida1 = cd_inducida(k1, CL)
    CD = CD01 + CD_inducida1  # Polar del avión.  Coeficiente de resistencia.
    '''
    -------------------------INICIO DE LA MANIOBRA-------------------------
    '''
    radius = v**2 / (g0 * (N - 1))  # Radio de giro (m).
    #Este radio de giro se obtiene para la velocidad inicial en vuelo
    # estacionario y para un factor de carga máximo según los pilones de carga
    # n = 3,5.
    dt = .1  # Diferencial de tiempo (s).
    '''------SISTEMA DE ECUACIONES PARA PRIMER TRAMO: VUELO ESTACIONARIO------
    '''
    #Energías.
    ecinetica = .5 * MASS * v**2  # Energía cinética (J).
    epotencial = MASS * g0 * h  # Energía potencial (J).
    emecanica = ecinetica + epotencial  # Energía mecánica (J).
    #Fuerzas.
    D = resistencia(v, rho, CD)  # Resistencia aerodinámica (N).
    L = sustentacion(v, rho, CL)  # Sustentación aerodinámica (N).
    Th = thrust(M, rho)  # Empuje (N).
    diferencia_T_D = Th - D
    #Esto nos va a permitir calcular en qué momento el empuje se verá superado
    # por la resistencia
    n = L / W  # Factor de carga.
    #Condiciones iniciales para la integración.
    t = 0
    x = 0
    omega = v / radius  # Velocidad angular en la maniobra de giro (rad/s).
    #Segunda ley de Newton en el eje horizontal (ejes cuerpo).
    dv = dt * (Th * cos(alfa) - D - MASS * g0 * sin(gama)) / MASS
    vx = v * cos(gama)  # Velocidad horizontal del avión en ejes tierra (m/s).
    vy = v * sin(gama)  # Velocidad vertical del avión en ejes tierra (m/s).
    dx = vx * dt  # Variación horizontal de la posición en ejes tierra (m).
    dh = vy * dt  # Variación vertical de la posición en ejes tierra (m).
    dtheta = v * dt / radius  # Variación del ángulo de empuje.
    vl = v
    Ml = vl / (GAMMA * R_AIR * T)**.5
    Cdl = cdll(Ml)
    D_misil = .5 * rho * Cdl * SREF_MISIL * vl**2
    '''-------SISTEMA DE ECUACIONES PARA SEGUNDO TRAMO: MANIOBRA DE GIRO-------
    Ahora comienza el bucle relativo al giro ascendente, que analiza la
    trayectoria con nuevas ecuaciones y condiciones de vuelo que se detallan
    más adlante.  El significado de theta < beta implica que el bucle realice
    el cálculo requerido siempre que el ángulo theta sea menor que beta.  Se ha
    obligado a que beta sea el ángulo de final de giro (al inicio del programa
    se le ha dado un rango de valores).  Por tanto, una vez que theta sea igual
    a beta, se dará por concluida la maniobra de giro y se comenzará con el
    siguiente tramo. Todo ello mientras la velocidad sea positiva, condición
    que permite ver si el resultado es realista o no.
    '''
    while gama < beta and v > 0:
        #A continuación, se escriben en un fichero todas las variables aquí
        # detalladas para cada valor de theta < beta.
        f.write('{0:.2f}\t'.format(t))  # Tiempo (s).
        f.write('{0:.3f}\t'.format(h))  # Altitud (m).
        f.write('{0:.3f}\t'.format(v))  # Velocidad (m/s).
        f.write('{0:.3f}\t'.format(M))  # Número de Mach.
        f.write('{0:.3f}\t'.format(alfa_grados))  # Ángulo de ataque (deg).
        f.write('{0:.3f}\t'.format(gama_grados))
        # Asiento de la velocidad (deg).
        f.write('{0:.3f}\t'.format(theta_grados))  # Ángulo de asiento (deg).
        f.write('{0:.3f}\t'.format(emecanica))  # Energía mecánica (J).
        f.write('{0:.3f}\t'.format(D))  # Fuerza de resistencia (N).
        #Ya que este análisis de maniobra, a diferencia del anterior, lleva un
        # cálculo para distintos valores de tiempo y velocidad, se debe
        # programar su evolución en términos de sus variaciones diferenciales.
        #Inicialización de variables y diferenciales para la maniobra del
        # misil.
        thetal = gama  # Inicialización del ángulo de asiento.
        thetalgrados = degrees(thetal)  # Ángulo de asiento (deg).
        yl = h  # Inicialización de la altitud (m).
        vl = v  # Inicialización de la velocidad (m).
        vxl = vl * cos(thetal)
        # Inicialización de la componente horizontal de velocidad (m/s).
        vyl = vl * sin(thetal)
        # Inicialización de la componente vertical de velocidad (m/s).
        tl = 0  # Inicialización temporal.
        xl = 0  # Inicialización de la posición en el eje x.
        sl = 0  # Inicialización del arco recorrido.
        dvxl = 0
        # Inicialización del diferencial de la componente horizontal de la
        # velocidad (m/s).
        dvyl = 0
        # Inicialización del diferencial de la componente vertical de la
        # velocidad (m/s).
        dsl = 0  # Inicialización del diferncial del arco recorrido (m).
        dxl = 0  # Inicialización del diferencial de la posición (m).
        dyl = 0  # Inicialización del diferencial de la altitud (m).
        dtl = .1  # Inicialización del diferencial de tiempo (m).
        dthetal = 0  # Inicialización del diferencial del ángulo de asiento.
        thetal = gama
        # El ángulo de asiento del avión es igual al ángulo de asiento de la
        # velocidad del misil.
        thetalgrados = degrees(thetal)
        Ddsl = 0
        #Aquí se escriben los incrementos diferenciales de las coordenadas
        # espaciales, temporales y de velocidad.  Los diferenciales se obtienen
        # del tramo anterior, y sirven para que estos nuevos valores actúen en
        # las nuevas condiciones para calcular nuevas variables.
        t = t + dt  # Evolución temporal (s).
        v = v + dv  # Velocidad (m/s).
        x = x + dx  # Posición horizontal (m).
        h = h + dh  # Altitud (m).
        r = RT + h  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        #Las variables termodinámicas habrán variado con la altura.
        rho = density(h)  # Densidad (kg/m3).
        T = temperature(h)  # Temperatura (K).
        Mu_Visc = viscosity(h)  # Viscosidad (Pa s).
        M = v / (GAMMA * R_AIR * T)**.5  # Mach de vuelo.
        n = 3.5  # Tomamos la condición de factor de carga máximo y constante.
        radius = v**2 / (g0 * (n - 1))  # Radio de giro varía con la velocidad.
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
        # Ángulos en grados para la exportación a los ficheros.
        alfa_grados = degrees(alfa)  # Ángulo de ataque (deg).
        theta_grados = degrees(theta)  # Ángulo de asiento (deg).
        gama_grados = degrees(gama)  # Asiento de la velocidad (deg).
        CL = alfa * CL_alfa1  # Coeficiente de sustentación.
        n = .5 * rho * v**2 * S_W * CL / W
        #Nuevas proyecciones de la velocidad (ejes tierra).
        vx = v * cos(gama)  # Proyección horizontal de la velocidad (m/s).
        vy = v * sin(gama)  # Proyección vertical de la velocidad (m/s).
        CD_inducida1 = cd_inducida(k1, CL)
        CD = CD01 + CD_inducida1  # Polar del avión.
        D = resistencia(v, rho, CD)  # Fuerza de resistencia (N).
        L = n * W  # Fuerza de sustentación (N).
        #Energías.
        ecinetica = .5 * MASS * v**2  # Energía cinética (J).
        epotencial = MASS * g0 * h  # Energía potencial (J).
        emecanica = ecinetica + epotencial  # Energía mecánica (J).
        #Empuje.
        Th = thrust(M, rho)  # Empuje (N).
        diferencia_T_D = Th - D
        #Velocidad angular.
        omega = v / radius  # Velocidad angular (rad/s).
        #Nuevas variaciones diferenciales.  Segunda ley de Newton.
        dv = dt * (Th * cos(alfa) - D - W * sin(gama)) / MASS
        #Esta última ecuación nos permite obtener un nuevo incvremento de
        # velocidad (m/s) que, al reiniciar el bucle será sumado al último
        # valor de la velocidad, obteniendo la variación de la velocidad.
        dx = vx * dt  # Variación horizontal de la posición (m).
        dh = v * sin(gama) * dt  # Variación vertical de la posición (m).
        dtheta = omega * dt  # Variación del ángulo de asiento.
        masa_misil = 1000  # Masa del misil (kg).
        #----------------PUESTA EN ÓRBITA DEL MISIL----------------
        while thetal > 0 and yl < 500000:
            #Con este bucle se calculan todas las variables  correspondientes
            # al lanzamiento del cohete.  La condición de parada del bucle es
            # que el ángulo de asiento del misil deje de ser positivo. Esto es
            # cuando el misil se encuentre completamente en posición horizontal
            # y que se alcancen los 500km de altitud.
            tl = tl + dtl  # Evolución temporal (s).
            xl = xl + dxl  # Posición horizontal (m).
            yl = yl + dyl  # Altitud (m).
            sl = sl + dsl  # Arco recorrido (m).
            thetal = thetal + dthetal  # Ángulo de asiento.
            thetalgrados = degrees(thetal)  # Ángulo de asiento (deg).
            vxl = vxl + dvxl  # Componente horizontal de la velocidad (m/s).
            vyl = vyl + dvyl  # Componente vertical de la velocidad (m/s).
            vl = (vxl**2 + vyl**2)**.5  # Módulo de la velocidad (m/s).
            r = RT + yl  # Distancia al centro de la Tierra (m).
            g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
            rho = density(yl)  # Densidad (kg/m3).
            T = temperature(yl)  # Temperatura (K).
            Mu_Visc = viscosity(yl)  # Viscosidad (Pa s).
            Ml = vl / (GAMMA * R_AIR * T)**.5  # Mach de vuelo.
            if tl <= t_combustion:
                RATIO_AREAS = 1 - SGASES / SREF_MISIL
            else:
                RARIO_AREAS = 1
            Cdl = cdll(Ml)  # Coeficiente de resistencia.
            D_misil = .5 * rho * Cdl * SREF_MISIL * vl**2
            # Fuerza de resistencia (N).
            Dx = D_misil * cos(thetal)
            # Componente horizontal de la fuerza de resistencia (N).
            Dy = D_misil * sin(thetal)
            # Componente vertical de la fuerza de resistencia (N).
            dvxl = -Dx / masa_misil * dtl
            # Diferencial de la componente horizontal de la velocidad (m/s).
            dvyl = -g0 * dtl - Dy / masa_misil * dtl
            # Diferencial de la componente vertical de la velocidad (m/s).
            if tl <= t_combustion:
                dvxl = dvxl + empuje_misil * cos(thetal) * dtl / masa_misil
                dvyl = dvyl + empuje_misil * sin(thetal) * dtl / masa_misil
                masa_misil = masa_misil - gasto * dtl
                # Variación de masa lineal (kg).
            dthetal = -dtl * g0 * cos(thetal) / vl
            # Diferencial del ángulo de asiento.
            dxl = vxl * dtl  # Diferencial de la posición horizontal (m).
            dyl = vyl * dtl  # Diferencial de la altitud (m).
            dsl = vl * dtl  # Diferencial del arco recorrido (m).
        Emec_misil = masa_misil * (GRAV * yl + vl**2 / 2)
        # Energía mecánica final del misil (J).
        # Trabajo de la resistencia (J).
        f.write('{0:.3f}\t'.format(tl))  # Tiempo (s).
        f.write('{0:.3f}\t'.format(yl))  # Altitud (m).
        f.write('{0:.3f}\t'.format(vl))  # Velocidad (m/s).
        f.write('{0:.3f}\t'.format(Ml))  # Mach de vuelo.
        f.write('{0:.3f}\t'.format(thetalgrados))  # Ángulo de asiento (deg).
        f.write('{0:.3f}\t'.format(masa_misil))  # Masa del misil (kg).
        f.write('{0:.3f}\t'.format(Emec_misil))  # Energía mecánica final (J).
        f.write('{0:.3f}\t'.format(xl))  # Posición x final (m).
        f.write('{0:.3f}\t'.format(empuje_misil))  # Empuje (N).
        f.write('{0:.3f}\n'.format(isp))  # Impulso específico (N s/kg).
    f.close()


#Como resumen:
# 1) El código ha empezado en una condición de vuelo uniforme.
# 2) La siguiente maniobra es un giro ascendente, a factor de carga máximo y
# constante, y con mínima resistencia.
# 3) Despliegue del cohete REOS. 1 Etapa. Con empuje y resistencia
# aerodinámica.
#
# Este programa exportará un archivo que, exportado a Excel nos permite
# observar cómo cambian las variables según las condiciones de vuelo.
