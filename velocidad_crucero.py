# -*- coding: utf-8 -*-
"""

@authrustor: Team REOS

Esta función proporciona la velocidad exacta de crucero para que el
empuje sea igual a la resistencia.

"""

from inputs_iniciales import Z0, MASS
from modelo_msise00 import GAMMA, R_AIR, temperature, density
from gravedad import gravity
from modelo_empuje import thrust
from aero_avion import k, cd0, cd_inducida, S_W, CD_interferencia
from aero_avion import resistencia


def vuelo_crucero(mach):

    '''Función que calcula la velocidad de crucero del avión. Los datos
    iniciales los toma del módulo "datos".
    '''

    temperature1 = temperature(Z0)  # Temperatura inicial.
    rho = density(Z0)  # Densidad inicial del aire (kg/m3).
    g_0 = gravity(Z0)  # Aceleración gravitatoria inicial (m/s2).
    weight = MASS * g_0
    # Peso inicial del avión dependiente de la gravedad
    # A la altura inicial el avión vuela en vuelo estacionario.
    velocity = mach * (GAMMA * R_AIR * temperature1)**.5
    # velocidad inicial (m/s).

    # Coeficientes aerodinámicos
    # Ecuacion de vuelo en crucero L = W
    c_l = 2 * weight / (rho * velocity**2 * S_W)
    k_1 = k(mach)
    c_d = cd0(mach) + cd_inducida(k_1, c_l)  # Polar del avión.

    cd_avion = c_d + CD_interferencia(mach)
    # Coeficiente de resistencia total.

    # Fuerzas.
    drag_aircraft = resistencia(velocity, rho, cd_avion)
    # Resistencia aerodinámica (N).
    thrust_aircraft = thrust(mach, rho)  # Empuje (N).
    # Esta es la ecuacion en eje horizontal T = D que es la condicion
    # que queremos cumplir por ello calculamos la diferencia y en el
    # while se intenta que sea 0.
    diferencia_t_d = abs(thrust_aircraft - drag_aircraft)
    # Se hace el absoluto para que no haya problemas.
    dmach = 0.000001  # Variaremos el Mach para iterar.

    while diferencia_t_d >= 10:

        mach = mach + dmach
        velocity = mach * (GAMMA * R_AIR * temperature1)**.5
        # velocidad inicial (m/s).

        # Coeficientes aerodinámicos
        # Ecuacion de vuelo en crucero L = W
        c_l = 2 * weight / (rho * velocity**2 * S_W)
        k_1 = k(mach)
        c_d = cd0(mach) + cd_inducida(k_1, c_l)  # Polar del avión.

        cd_avion = c_d + CD_interferencia(mach)
        # Coeficiente de resistencia.

        # Fuerzas.
        drag_aircraft = resistencia(velocity, rho, cd_avion)
        # Resistencia aerodinámica (N).
        thrust_aircraft = thrust(mach, rho)  # Empuje (N).
        # Esta es la ecuacion en eje horizontal T = D que es la
        # condicion que queremos cumplir por ello calculamos la
        # diferencia y en el While se intenta que sea 0.

        diferencia_t_d = abs(thrust_aircraft - drag_aircraft)

    print('El número de mach es: ', mach)

    return mach
