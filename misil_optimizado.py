# -*- coding: utf-8 -*-
"""

@author: Team REOS

Programa que nos permite calcular las masas del cohet a lanzar
desde los distintos puntos de lanzamiento, optimizandolo para obtener
la mínima masa total del cohete para una carga de pago dada o viceversa.
Los cohetes optimizados son de 2 etapas y se lanzan desde los puntos
obtenidos por el modulo avion,

"""
from math import radians, cos, sin, degrees, pi, exp

from modeloISA import density, temperature, GAMMA, viscosity, R_AIR, gravity
from aero_misil import cdll, SREF_MISIL, SGASES
from avion import velocidad, altura, psi_list, fi_list

# --------------------------CONDICIONES GRAVITATORIAS-------------------
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2).
FI_GRADOS = fi_list
PSI_GRADOS = psi_list

MASA_PAYLOAD = 34
for gasto in [14]:
    # Creamos un string a partir del gasto para crear archivos con el
    # nombre del valor de la variable gasto.
    gastotexto = str(gasto)
    # Generación de 'archivo', donde se almacenaran los datos.
    # Cada etiqueta verde es el titulo de cada columna de valores.
    archivo = open(gastotexto, "w")
    archivo.write('Time \tAltura lanzamiento \tAngulo lanz'
                  '\tVelocidad final \tVelocidad perdida'
                  '\tmasa total \tfil  \tindice \n')
    #
    # Este bucle permite recorrer los elementos de los vectores creados
    # en el archivo avión, velocidad, altura, fi y psi.
    # Esto quiere decir que cada de los vectores
    # velocidad, altura, FI_GRADOS y PSI_GRADOS, indican
    # el un punto de lanzamiento calculado y asi con cada
    # indice, de tal manera que cada uno represneta un punto
    # de lanzamiento.
    for i in range(144, 284, 1):

        gasto_etapa2 = gasto
        isp1 = 280
        Empuje_misil = gasto * isp1 * 9.81
        # Relacion estrutural un 5% del peso del propulsante de etapa.
        ratio_estructural = .05
        # Incremento de un 6% de la masa final
        # para la maniobras de circulacion.
        masa_util = MASA_PAYLOAD * 1.06
        # Velocidades implicadas y parametros del optimizador
        v_orb = 7800
        v_loss = 782        # Velocidad con la que empieza a iterar
        v_loss_empirica = 0
        v_rotacional = 400
        error = 800
        contador = 1
        #
        # ----Iterador----
        #
        # Con este bucle realizamos la interación de la pérdida por
        # velocidad. A partir de un valor incial dado calcula la
        # trayectoria del misi y con ella obtenemos la verdadera
        # pérdida por velocidad. Tras ello se comparan valores.
        # Este proceso durará mientras que la diferencia entre ambos
        # valores sea mayor a 10, este valor se puede modificar para
        # mayor precisión, y mientras el contador de iteraciones sea
        # menor de 50, esto sirve para que si un caso no tiene solución
        # u oscila entre 2 valores mayores al error sirva de parada.
        while error > 10 and contador <= 50:

            # Los valores siguientes los obtiene de vectores generados
            # en el archivo avión con todos los puntos de lanzamiento.
            # Esto quiere decir que el indice 1 de los vectores
            # velocidad, altura, FI_GRADOS y PSI_GRADOS, indican
            # el primer punto de lanzamiento calculado y asi con cada
            # indice, de tal manera que cada uno represneta un punto
            # de lanzamiento.

            z0 = altura[i]  # Altura a la que se lanza.
            v0 = velocidad[i]  # Velocidad inicial con la que se lanza.
            fil_grados0 = FI_GRADOS[i]
            # Ángulo de ataque de la velocidad sobre la vertical local.
            psil_grados0 = PSI_GRADOS[i]
            # Ángulo de la vertical local en ejes centro tierra.
            theta_grados0 = 90 - fil_grados0 - psil_grados0

            # Asignamos a 'zl' el valor z0, el inicial, y será 'zl'
            # la varibale que varie durante la maniobra, consevando asi
            # el valor de lanzamiento inicial en 'z0'.
            zl = z0

            # Lo mismo con los ángulos que ocn 'zl', creamos variables
            # para no perder las variables de las condicones iniciales.
            fil_grados = fil_grados0
            psil_grados = psil_grados0
            fil = radians(fil_grados)  # Conversión a grados del ángulo.
            psil = radians(psil_grados)  # Conversión a grados de psil.

            rl = RT + zl  # Distancia al centro de la Tierra.

            # De nuevo creamos 'vl'. que ira variando a lo largo
            # de la maniobra y además añadimos el efecto de la velocidad
            # rotacional.
            vl = v0 + v_rotacional
            Ml = vl / ((GAMMA * R_AIR * temperature(zl))**.5)
            vxl = vl * sin(fil)
            # Inicialización de la componente horizontal de velocidad.
            vyl = vl * cos(fil)
            # Inicialización de la componente verical de velocidad.
            tl = 0
            xl = 0
            # Iniciamos la variable x de desplazamiento del misil como
            # el que lleva el avion en el punto de lanzamiento.
            sl = 0  # Inicialización del arco recorrido.
            dvl = 0
            dsl = 0  # Inicialización del diferncial del arco recorrido.
            dxl = 0  # Inicialización del diferencial de la posición.
            dzl = 0
            # Inicialización del diferencial de la altura del misil.
            dfil = 0
            # Inicialización de la variación del ángulo fi inicial.
            dpsil = 0
            # Inicialización de la variación del ángulo de psi inicial.
            dtl = 0.01
            # Diferencial de tiempo empleado en la integración.

            Cdl = cdll(Ml, zl)  # Coeficiente resistencia del misil
            D_misil = .5 * density(zl) * Cdl * SREF_MISIL * vl**2
            Emec_misil = 1000 * (gravity(zl) * zl + (vl ** 2) / 2)
            g = gravity(zl)  # Aceleración de la gravedad (m/s2).

            # ---------OPTIMIZADOR---------
            # Este optimizador de 2 etapas, nos permite calcular el
            # cohete de menor masa para una velocidad de orbita dada.
            # La otra incógnita sería la pérdida por velcoidad que
            # se calcula posteriomente mediante la interación en la
            # que se encuentra sumergido el optimizador.
            #
            # ----Parámetros implicados----
            #
            # Velocidad ideal, es la velocidad ideal que habría que
            # incrementar v0 para obtener la velcoidad orbital sin tener
            # en cuenta las perdidas por velocidad ni la velocidad
            # rotacional, pero no es el incremento que hayq ue dar.
            # Hay que tener en cuenta los demás efectos.

            velocidad_ideal = v_orb - v_rotacional + v_loss - v0
            # Velocidad ideal, necesaria para 'f'.
            f = (((1 + ratio_estructural) / exp(velocidad_ideal /
                  (2 * isp1 * gravity(zl)))) - ratio_estructural)
            # Coeficiente dado por el optimizador.

            # MASAS DEL COHETE

            masa_total = masa_util / (f**2)  # Masa inicial del cohete.
            masa_etapa2 = masa_util / f  # Masa de la segunda etapa.
            masa_propulsante_etapa1 = ((masa_total - masa_etapa2) /
                                       (1 + ratio_estructural))
            masa_propulsante_etapa2 = ((masa_etapa2 - masa_util) /
                                       (1 + ratio_estructural))
            masa_estructura1 = masa_propulsante_etapa1 * ratio_estructural
            masa_estructura2 = masa_propulsante_etapa2 * ratio_estructural

            t_combustion = masa_propulsante_etapa1 / gasto
            # Tiempo de combustión primera etapa.
            t_combustion_2_etapa = masa_propulsante_etapa2 / gasto_etapa2
            # Tiempo de combustión segunda etapa.
            retardo_encendido = 0
            # Tiempo desde que se apaga la primera etapa y se enciende
            # la segunda.
            t_fin_combustion2 = (t_combustion + t_combustion_2_etapa +
                                 retardo_encendido)
            # Instante temporal en el que se apaga la segunda etapa.
            masa_misil = masa_total
            # Esta masa irá variando a lo largo de la maniobra.
            v_loss_empirica = 0
            # Reinicio del valor de la pérdida de velocidad calculada,
            # asi al iniciar una nueva iteración vuelve a adoptar
            # el valor de 0 en el momento de iniciar la maniobra.

            # Este bucle indica la combustión de la primera etapa,
            # finalizará al llegar al tiempo de comuistión o si
            # se alcanza la condicion de que el misil cae.

            while tl <= t_combustion and fil < pi / 2:

                tl = tl + dtl  # Evolución temporal (s.
                xl = xl + dxl  # Posición horizontal (m).
                zl = zl + dzl  # Altitud (m).
                sl = sl + dsl  # Arco recorrido (m).
                rl = RT + zl  # Altura al centro de la Tierra.

                fil = fil + dfil
                # Variación del ángulo de asiento en vertical local.
                fil_grados = degrees(fil)
                # Conversión a grados de fil.
                psil = psil + dpsil
                # Variación del ángulo de la vertical local.
                psil_grados = degrees(psil)
                # Conversión a grados de psil.

                vl = vl + dvl  # Variación de la velocidad del misil.
                vxl = vl * sin(fil)  # Componente horizontal.
                vyl = vl * cos(fil)  # Componente vertical.

                # Variables terodinamicas.
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl)  # Viscosidad.
                gl = gravity(zl)  # Gravedad variable con altura.
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5)

                Ratio_areas = (SREF_MISIL - SGASES) / SREF_MISIL
                # Término que determina la resistencia de base,
                # al haber combustión deja de valer 1.
                Cdl = cdll(Ml, zl)
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2

                # Ecuación horizontal en ejes velocidad.
                # Ecuación pricipal que nos determina la velocidad,
                # a partir de ella se determinan el resto de variables.
                dvl = dtl / masa_misil * (Empuje_misil - D_misil
                                          - gl * masa_misil * cos(fil))

                masa_misil = masa_misil - gasto * dtl
                # Variación de masa en cada 'dt' debido a la combustión,
                dfil = gl / vl * dtl
                # Diferencial del ángulo de asiento.
                dpsil = vxl / rl * dtl
                # Diferencial de la vertical local.

                dv_loss = dtl / masa_misil * (D_misil + gl * masa_misil
                                              * cos(fil))
                # Diferencial de pérdidas por velocidad.
                v_loss_empirica = v_loss_empirica + dv_loss
                # Sumatorio de los diferenciales de pérdidas por
                # velocidad dan lugar a 'v_loss_empírica' que será
                # la velocidad que compararemos con la 'v_loss' iniciada
                # para iterar.

                dxl = vxl * dtl  # Diferencial de la posición (m).
                dzl = vyl * dtl  # Diferencial de la altitud (m).
                dsl = vl * dtl  # Diferencial del arco recorrido (m).

            masa_misil = masa_misil - masa_estructura1
            # Decremento de masa al soltarse la primera etapa.

            # Este bucle indica la combustión de la segunda etapa,
            # finalizará al llegar al tiempo de comuistión o si se
            # alcanza la condicion de que el misil cae o llega a 500km.

            while t_fin_combustion2 >= tl > t_combustion and fil < pi / 2 and zl < 500000:

                tl = tl + dtl  # Evolución temporal (s).
                xl = xl + dxl  # Posición horizontal (m).
                zl = zl + dzl  # Altitud (m) .
                sl = sl + dsl  # Arco recorrido (m).
                rl = RT + zl   # Altura al centro de la Tierra (m).

                fil = fil + dfil
                # Variación del ángulo de asiento en vertical local.
                fil_grados = degrees(fil)
                # Conversión a grados de fil.
                psil = psil + dpsil
                # Variación del ángulo de la vertical local.
                psil_grados = degrees(psil)
                # Conversión a grados de psil.

                vl = vl + dvl  # Variación de la velocidad del misil.
                vxl = vl * sin(fil)  # Componente horizontal.
                vyl = vl * cos(fil)  # Componente vertical.

                # Variables terodinamicas.
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl)  # Viscosidad.
                gl = gravity(zl)  # Gravedad variable con altura.
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5)  # Mach de vuelo,

                Ratio_areas = (SREF_MISIL - SGASES) / SREF_MISIL
                # Término que determina la resistencia de base,
                # al haber combustión deja de valer 1.
                Cdl = cdll(Ml, zl)
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2

                # Ecuación horizontal en ejes velocidad.
                # Ecuación pricipal que nos determina la velocidad,
                # a partir de ella se determinan el resto de variables.
                dvl = dtl / masa_misil * (Empuje_misil - D_misil
                                          - gl * masa_misil * cos(fil))

                masa_misil = masa_misil-gasto_etapa2 * dtl
                # Variación de masa en cada 'dt' debido a la combustión,
                dfil = gl / vl * dtl
                # Diferencial del ángulo de asiento.
                dpsil = vxl / rl * dtl
                # Diferencial de la vertical local.

                dv_loss = dtl / masa_misil * (D_misil + gl * masa_misil
                                              * cos(fil))
                # Diferencial de pérdidas por velocidad.
                v_loss_empirica = v_loss_empirica + dv_loss
                # Sumatorio de los diferenciales de pérdidas por
                # velocidad dan lugar a 'v_loss_empírica' que será
                # la velocidad que compararemos con la 'v_loss' iniciada
                # para iterar.

                dxl = vxl * dtl  # Diferencial de la posición (m).
                dzl = vyl * dtl  # Diferencial de la altitud (m).
                dsl = vl * dtl  # Diferencial del arco recorrido (m).

            masa_misil = masa_misil - masa_estructura2
            # Decremento de masa al soltarse la segunda etapa.

            while fil < pi / 2 and zl < 500000:

                tl = tl + dtl  # Evolución temporal (s).
                xl = xl + dxl  # Posición horizontal (m).
                zl = zl + dzl  # Altitud (m) .
                sl = sl + dsl  # Arco recorrido (m).
                rl = RT + zl   # Altura al centro de la Tierra (m).

                fil = fil + dfil
                # Variación del ángulo de asiento en vertical local.
                fil_grados = degrees(fil)
                # Conversión a grados de fil.
                psil = psil + dpsil
                # Variación del ángulo de la vertical local.
                psil_grados = degrees(psil)
                # Conversión a grados de psil.

                vl = vl + dvl  # Variación de la velocidad del misil.
                vxl = vl * sin(fil)  # Componente horizontal.
                vyl = vl * cos(fil)  # Componente vertical.

                # Variables terodinamicas.
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl)  # Viscosidad.
                gl = gravity(zl)  # Gravedad variable con altura.
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5)  # Mach de vuelo,

                Ratio_areas = 1
                # Término que determina la resistencia de base,
                # al no haber combustión deja de vale 1.
                Cdl = cdll(Ml, zl)
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2

                # Ecuación horizontal en ejes velocidad.
                # Ecuación pricipal que nos determina la velocidad,
                # a partir de ella se determinan el resto de variables.
                dvl = dtl / masa_misil * (- D_misil - gl * masa_misil
                                          * cos(fil))

                dfil = gl / vl * dtl
                # Diferencial del ángulo de asiento.
                dpsil = vxl / rl * dtl
                # Diferencial de la vertical local.

                dv_loss = - dvl
                # Diferencial de pérdidas por velocidad. En este caso
                # coincide con el diferencial de 'dvl', ya que solo se
                # pierde velocidad al no haber empuje.
                v_loss_empirica = v_loss_empirica + dv_loss
                # Sumatorio de los diferenciales de pérdidas por
                # velocidad dan lugar a 'v_loss_empírica' que será
                # la velocidad que compararemos con la 'v_loss' iniciada
                # para iterar.

                dxl = vxl * dtl  # Diferencial de la posición (m).
                dzl = vyl * dtl  # Diferencial de la altitud (m).
                dsl = vl * dtl  # Diferencial del arco recorrido (m).

            error = abs(v_loss_empirica - v_loss)
            v_loss = v_loss_empirica

            if contador == 50:

                print('error de interacion en el ángulo ', fil_grados0,
                      'con un error de ', error, ' m/s')
            contador = contador + 1

        print('para ', theta_grados0, 'la masa del misil es ', masa_total,
              ' la altura es ', zl, ' y el angulo de llegada', fil_grados)
        print('El gasto es: ', gasto, 'el numeor de la lista es ', i)

        archivo.write('%.8f\t' % gasto)
        archivo.write('%.8f\t' % z0)
        archivo.write('%.8f\t' % theta_grados0)
        archivo.write('%.8f\t' % zl)
        archivo.write('%.8f\t' % vl)
        archivo.write('%.8f\t' % v_loss)
        archivo.write('%.8f\t' % masa_total)
        archivo.write('%.8f\t' % fil_grados)
        archivo.write('%.8f\n' % i)

        if zl > 500000:
            break

    archivo.close()
