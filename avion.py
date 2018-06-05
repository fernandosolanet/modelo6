# -*- coding: utf-8 -*-
"""
@author: Team REOS
"""

from math import radians, cos, sin, degrees, pi

import sys
sys.path.insert(0, '/path/to/modulos')
sys.path.insert(0, '/path/to/aerodinamica')
sys.path.insert(0, '/path/to/atmosfera')
sys.path.insert(0, '/path/to/empuje')


from inputs_iniciales import MASS
from modulos.atmosfera.modelo_msise00 import GAMMA, viscosity, R_AIR
from modulos.atmosfera.modelo_msise00 import density, temperature, pressure
from modulos.atmosfera.gravedad import gravity, RT

from modulos.empuje.modelo_empuje import thrust
from modulos.aerodinamica.aero_avion import cl, angulo_ataque, angulo_ataqueT, S_W, k, cd_inducida
from modulos.aerodinamica.aero_avion import cd_interferencia, cd0
from modulos.aerodinamica.aero_avion import resistencia, sustentacion






lanz = open('lanzamientos.txt', 'r')

#lanz = open('lanzamientos.txt', 'r')

altura_inicio  = []
mach_inicio = []

for line in lanz:
    s = line.strip()
    s = s.split()
    altura_inicio.append(float(s[0]))  # Altura .
    mach_inicio.append(float(s[1]))  # Mach
lanz.close()

DT = 0.1  # Diferencial de tiempo (s).

a = 0
b = len(mach_inicio)

for i in range(a, b):
  if __name__ == '__main__':
      N = 3.5  # Factor de carga máximo.
      BETA = radians(89)
      altura = str(round(altura_inicio[i],0))
      nombre = str(altura+' m')
      F = open(nombre, 'w')  # Fichero de escritura sin extensión.

      # Las variables que se escriben a continuación corresponden al paso de
      # integración del avión en cada paso.

      F.write('TIME (s)\tALTURA lanzamiento (m) \tDesplazamiento (m)')
      F.write('\tVELOCIDAD (m/s)\tMACH\tD (N)\tALFA (deg)\tGAMMA (deg)')
      F.write('\tTHETA (deg) \tfi (deg) \tpsi (deg)\n')

      # Las variables que es escriben en el archivo que aparecen a
      # continuación se corresponden con el punto final de maniobra del misil
      # para cada lanzamiento en cada paso de integración del misil

      # Cabezas de tabla.

      # -------------------------CONDICIONES INICIALES-----------------------
      # Ahora, para los próximos cálculos, se definen las variables
      # termodinámicas obtenidas del modelo atmosférico.

      VELOCIDAD = []
      ALTURA = []
      FI_LIST = []
      PSI_LIST = []

      # Vectores VELOCIDAD, ALTURA, fil y psil vacios, en ellos guardaremos
      # los valores de las distintas variables para su posterior uso

      M = mach_inicio[i]
      Z0 = altura_inicio[i]
      R = RT + Z0  # Distancia desde el centro de la Tierra (m)
      G0 = gravity(Z0)  # Aceleración gravitatoria inicial (m/s2)
      RHO = density(Z0)  # Densidad inicial del aire (kg/m3)
      P = pressure(Z0)  # Presión inicial del aire (Pa)
      TEMPERATURE = temperature(Z0)  # Temperatura inicial del aire (K)
      MU_VICS = viscosity(Z0)  # Viscosidad inicial
      W = MASS * G0  # Peso inicial del avión dependiente de la gravedad

      # A la ALTURA inicial el avión vuela en vuelo estacionario
      
      V = M * (GAMMA * R_AIR * TEMPERATURE)**.5  # Velocidad inicial (m/s)

      CL1 = 2 * W / (RHO * V**2 * S_W)
      K1 = k(M)
      CD01 = cd0(M)
      CD_INDUCIDA1 = cd_inducida(K1, CL1)
      CD_INTERFERENCIA1 = cd_interferencia(M)
      CD = CD01 + CD_INDUCIDA1  # Polar del avión.  Coeficiente de resistencia
      CD_AVION = CD + 0.35*CD_INTERFERENCIA1
      
      # Ángulos de asiento, de ataque y de asiento de la VELOCIDAD iniciales.
      ALFA_NUMERICO = angulo_ataque(M, CL1)
      ALFA_NUMERICO_GRADOS = degrees(ALFA_NUMERICO)
      ALFA = angulo_ataqueT(ALFA_NUMERICO, M)  # Ángulo de ataque.
      ALFA_GRADOS = degrees(ALFA)  # Ángulo de ataque en grados.
      GAMA = 0  # Ángulo de asiento de la VELOCIDAD.
      GAMA_GRADOS = degrees(GAMA)
      # Ángulo de asiento de la VELOCIDAD en grados.
      THETA = GAMA + ALFA  # Ángulo de asiento de la VELOCIDAD.
      THETA_GRADOS = degrees(THETA)  # Ángulo de asiento en grados.

      # Recálculo de CL con el ángulo de ataque.
      CL1 = cl (M, ALFA)
      CD_INDUCIDA1 = cd_inducida(K1, CL1)
      CD = CD01 + CD_INDUCIDA1  # Polar del avión.  Coeficiente de resistencia.
      CD_AVION = CD + 0.35*CD_INTERFERENCIA1

      # Fuerzas
      D = resistencia(V, RHO, CD_AVION)  # Resistencia aerodinámica (N)

      TH = thrust(M, RHO)  # Empuje (N)
      L = sustentacion(V, RHO, CL1)

      # Esta es la ecuación en eje horizontal T = D que es la condicion que
      # queremos cumplir por ello calculamos la diferencia y en el while se
      # intenta que sea 0

      DIFERENCIA_T_D = TH - D

      # Ángulos de movimiento bidimensional terrestre iniciales
      PSI = 0  # Ángulo desplazado de la Tierra
      FI = - PSI - GAMA + pi / 2
      # Ángulo de la VELOCIDAD sobre vertical local
      PSI_GRADOS = degrees(PSI)  # Ángulo desplazado de la Tierra en grados
      FI_GRADOS = degrees(FI)
      # Ángulo velocidad sobre vertical local en grados

      VELOCIDAD.append(V)
      PSI_LIST.append(PSI_GRADOS)
      ALTURA.append(Z0)
      FI_LIST.append(FI_GRADOS)
      # ------------------------INICIO DE LA MANIOBRA------------------------

      # Se impone un valor constante del radio de giro, es decir, giro
      # ascendente a radio constante.  En futuras versiones de este cálculo,
      # esto se cambiará para buscar una solución más realista.  Por ahora,
      # con objeto de facilitar cálculos, se dejará así.

      RADIUS = V**2 / (G0 * (N - 1))  # Radio de giro inicial (m).
      # Este radio de giro se obtiene para la VELOCIDAD inicial en vuelo
      # estacionario y para un factor de carga máximo según los pilones de
      # carga, n = 3,5.

      # -----SISTEMA DE ECUACIONES PARA PRIMER TRAMO: VUELO ESTACIONARIO------

      # Energías
      E_CINETICA = .5 * MASS * V**2  # Energía cinética (J)
      E_POTENCIAL = MASS * G0 * Z0  # Energía potencial (J)
      E_MECANICA = E_CINETICA + E_POTENCIAL  # Energía mecánica (J)

      N = round(L / W,0)  # Factor de carga.

      # Condiciones iniciales para la integración
      T = 0
      X = 0
      Z = Z0
      XT = X
      ZT = Z0
      OMEGA = V / RADIUS  # Velocidad angular en la maniobra de giro (rad/s)

      # Segunda ley de Newton en el eje horizontal (ejes cuerpo)
      DV = DT * (TH * cos(ALFA) - D - W * cos(FI)) / MASS

      VX = V * sin(FI)  # Velocidad horizontal del avión en ejes HL
      VY = V * cos(FI)  # Velocidad vertical del avión en ejes VL

      VXT = V * cos(GAMA)  # Veloc horizontal del avión en ejes tierra
      VYT = V * sin(GAMA)  # Veloc vertical del avión en ejes tierra

      DX = VX * DT + DV * cos(FI) * DT / 2
      # Variación horizontal de la posición en "Horizonte local"
      DZ = VY * DT + DV * sin(FI) * DT / 2
      # Variación vertical de la posición en "Vertical local"

      DXT = VXT * DT  # Variación horizontal de la posición en ejes tierra
      DZT = VYT * DT  # Variación vertical de la posición en ejes tierra

      DTHETA = V * DT / RADIUS  # Variación del ángulo de asiento
      DPSI = V * sin(FI) / R * DT  # Variación del ángulo de la vertical local

      # ------SISTEMA DE ECUACIONES PARA SEGUNDO TRAMO: MANIOBRA DE GIRO------
      # Ahora comienza el bucle relativo al giro ascendente, que analiza la
      # trayectoria con nuevas ecuaciones y condiciones de vuelo que se
      # detallan más afelante.  El significado de THETA < BETA implica que el
      # bucle realice el cálculo requerido siempre que el ángulo THETA sea
      # menor que BETA.  Se ha obligado a que BETA sea el ángulo de final de
      # giro (al inicio del programa se le ha dado un rango de valores).
      # Por tanto, una vez que THETA sea igual a BETA, se dará por concluida
      # la maniobra de giro y se comenzará con el siguiente tramo. Todo ello
      # mientras la velocidad sea positiva, condición que permite ver si el
      # resultado es realista o no.

      while GAMA < BETA and V > 0:
          # A continuación, se escriben en un fichero todas las variables aquí
          # detalladas para cada valor de THETA < BETA.

          # Ya que este análisis de maniobra, a diferencia del anterior
          # lleva un cálculo para distintos valores de tiempo y velocidad,
          # se debe programar su evolución en términos de sus variaciones
          # diferenciales.  Aquí se escriben los incrementos diferenciales
          # de las coordenadas espaciales, temporales y de velocidad.
          # Los diferenciales se obtienen del tramo anterior, y sirven para
          # que estos nuevos valores actúen en las nuevas condiciones para
          # calcular nuevas variables.

          F.write('%.8f\t' % T)  # Tiempo (s)
          F.write('%.8f\t' % Z)  # Altitud (m)
          F.write('%.8f\t' % X)  # Recorrido (m)
          F.write('%.8f\t' % V)  # Velocidad (m/s)
          F.write('%.8f\t' % M)  # Número de Mach
          F.write('%.8f\t' % D)  # Resistencia
          F.write('%.8f\t' % ALFA_GRADOS)  # Ángulo de ataque (grados)
          F.write('%.8f\t' % GAMA_GRADOS)  # Asiento de la VELOCIDAD (grados)
          F.write('%.8f\t' % THETA_GRADOS)  # Ángulo de asiento (grados)
          F.write('%.8f\t' % FI_GRADOS)     # Ángulo del avión sobre la VL
          F.write('%.8f\t' % PSI_GRADOS)    # Ángulo de la vertical local
          F.write('%.8f\n' % ALFA_NUMERICO)  # Ángulo de ataque (grados)

          # Inicio del cálculo de la maniobra del avión

          # Incrementos
          T = T + DT  # Evolución temporal (s).
          V = V + DV  # Velocidad (m/s).
          X = X + DX  # Posición horizontal-local (m).
          XT = XT + DXT  # Posición horizontal en ejes tierra (m)
          Z = Z + DZ  # Altitud en horizonte-local (m).
          ZT = ZT + DZT  # Altitud en ejes del punto de lanzamiento.

          # Guardamos las variables en la lista
          VELOCIDAD.append(V)
          ALTURA.append(Z)

          # Gravitación
          R = RT + Z  # Distancia al centro de la Tierra (m).
          G = gravity(Z)  # Aceleración de la gravedad (m/s2).
          W = MASS * G  # Peso dependiente de la gravedad.

          # Las variables termodinámicas habrán variado con la ALTURA.
          RHO = density(Z)  # Densidad (kg/m3).
          TEMPERATURE = temperature(Z)  # Temperatura (K).
          MU_VICS = viscosity(Z)  # Viscosidad

          M = V / (GAMMA * R_AIR * TEMPERATURE)**.5  # Mach de vuelo.
          N = 3.5  # Tomamos la condición de factor de carga máximo y cte.
          RADIUS = V**2 / (G * (N - 1))  # Radio de giro varía con la veloc.

          # Las características aerodinámicas varían con el número de Mach.
          # Ecuación eje vertical.
          CL1 = 2 * W *N / (RHO * V**2 * S_W)
          
          # Como las gráficas no tienen valores mayores de 1, hay que imponer
          # un ángulo de ataque que asegure que el avión en pérdida, es decir 
          # CL > 1 implica ángulos de ataque muy altos que no están implementados
          #por falta de información. 
          if CL1 > 1.01691651:
              ALFA_NUMERICO = radians(50)
          else:
              ALFA_NUMERICO = round(angulo_ataque(M, CL1),5)
              
          # El nuevo ángulo de ataque resultará del nuevo factor de carga
          # (antes 1 y ahora máximo), la nueva VELOCIDAD y las nuevas
          # características aerodinámicas.
#          ALFA_NUMERICO = round(angulo_ataque(M, CL1),5)
          ALFA_NUMERICO_GRADOS = degrees(ALFA_NUMERICO)
          ALFA = round(angulo_ataqueT(ALFA_NUMERICO, M),5)  # Ángulo de ataque.
          
          K1 = k(M)
          CD01 = cd0(M)
          CD_INDUCIDA1 = cd_inducida(K1, CL1)
          CD = CD01 + CD_INDUCIDA1  # Polar del avión.  Coeficiente de resistencia.
          CD_INTERFERENCIA1 = cd_interferencia(M)
          CD_AVION = CD + 0.35*CD_INTERFERENCIA1

          THETA = THETA + DTHETA  # Ángulo de asiento (empuje horizontal)
          GAMA = THETA - ALFA  # Ángulo de asiento de la velocidad
          PSI = PSI + DPSI  # Ángulo de vertical local
          FI = - PSI - GAMA + pi / 2
          # Ángulo de la vel sobre la vertical local

          # Ángulos en grados para la exportación a los ficheros
          ALFA_GRADOS = degrees(ALFA)  # Ángulo de ataque (grados)
          THETA_GRADOS = degrees(THETA)  # Ángulo de asiento (grados)
          GAMA_GRADOS = degrees(GAMA)  # Asiento de la VELOCIDAD (grados)
          PSI_GRADOS = degrees(PSI)    # Ángulo de vertical local (grados)
          FI_GRADOS = degrees(FI)
          # Ángulo velocidad sobre la vertical local (grados)

          # Guardamos las variables en la lista
          PSI_LIST.append(PSI_GRADOS)
          FI_LIST.append(FI_GRADOS)

          #Se recalcula el CL en función de si el avión entra en pérdida o no
          # Se recalcula n con la ecuación del eje vertical
          if ALFA_NUMERICO != ALFA:

            K1 = k(M)
            if CL1 > 1.01691651:
                       CL1 = 1.01691651
            CD01 = cd0(M)
            CD_INDUCIDA1 = cd_inducida(K1, CL1)
            CD = CD01 + CD_INDUCIDA1  # Polar del avión.  Coeficiente de resistencia.
            CD_INTERFERENCIA1 = cd_interferencia(M)
            CD_AVION = CD + 0.35*CD_INTERFERENCIA1
            
            N = .5 * RHO * V**2 * S_W * CL1 / W

          # Nuevas proyecciones de la VELOCIDAD (ejes radiales y ejes
          # posición inicial)

          VX = V * sin(FI)  # Velocidad horizontal en ejes locales
          VY = V * cos(FI)  # Velocidad vertical en ejes locales
          VXT = V * cos(GAMA)  # Proyección horiz velocidad en ejes tierra
          VYT = V * sin(GAMA)  # Proyección vertical velocidad en ejes tierra
          
          D = resistencia(V, RHO, CD_AVION)  # Fuerza de resistencia (N)
          L = N * W  # Fuerza de sustentación (N)

          # Energías
          E_CINETICA = .5 * MASS * V**2  # Energía cinética (J)
          E_POTENCIAL = MASS * G * Z  # Energía potencial (J)
          E_MECANICA = E_CINETICA + E_POTENCIAL  # Energía mecánica (J)

          # Empuje
          TH = thrust(M, RHO)  # Empuje (N)
          DIFERENCIA_T_D = TH - D  # Diferencia de empuje frente a resitencia

          # Velocidad angular
          OMEGA = V / RADIUS  # Velocidad angular (rad/s)

          # Nuevas variaciones diferenciales.  Segunda ley de Newton.
          DV = DT * (TH * cos(ALFA) - D - W * cos(FI)) / MASS
          DPSI = ((V * sin(FI)) / (R)) * DT
          # Esta última ecuación nos permite obtener un nuevo incvremento de
          # velocidad (m/s) que, al reiniciar el bucle será sumado al último
          # valor de la velocidad, obteniendo la variación de la velocidad.

          # Diferenciales
          DX = VX * DT + DV * cos(FI) * DT / 2
          # Variación horizontal de la posición (m)
          DZ = VY * DT + DV * sin(FI) * DT / 2
          # Variación vertical de la posición (m)

          DXT = VXT * DT  # Variación horiz. posición en ejes tierra (m/s)
          DZT = VYT * DT  # Variación vertical posición en ejes tierra (m/s)

          DTHETA = OMEGA * DT  # Variación del ángulo de asiento

      F.close()

      E = open('exports', 'w')
      for i, a in enumerate(ALTURA):
          E.write(str(a) + ' ')

      E.write('\n')
      for i, a in enumerate(VELOCIDAD):
          E.write(str(a) + ' ')

      E.write('\n')
      for i, a in enumerate(FI_LIST):
          E.write(str(a) + ' ')

      E.write('\n')
      for i, a in enumerate(PSI_LIST):
          E.write(str(a) + ' ')

      E.close()

# Como resumen:
# 1) El código ha empezado en una condición de vuelo uniforme.
# 2) La siguiente maniobra es un giro ascendente, a factor de carga
# máximo y constante, y con mínima resistencia.
# 3) La última maniobra es un ascenso con el ángulo final del giro, con
# coeficiente de sustentación óptimo.
# Este programa nos exportará un total de 9 archivos (uno para cada
# incremento de 10º del ángulo final de giro) que, exportados a Excel
# nos permiten observar cómo cambian las variables según las condiciones
# de vuelo.
