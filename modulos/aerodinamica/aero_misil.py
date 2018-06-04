# -*- coding: utf-8 -*-
"""
@author: Team REOS
Funciones que nos dan los coeficientes aerodinámicos del avión F4.
Coeficientes aerodinámicos de sustentación y de resistencia y todo
aquello relacionado con los mismos.
"""

from math import log10, pi, degrees, atan, sqrt

from modulos.atmosfera.modelo_msise00 import GAMMA, density, temperature, R_AIR, viscosity
from inputs_iniciales import (DIAMETRO_M, LONGITUD_CONO, LONGITUD_MISIL,
                              NUM_ALETAS, ESPESOR_ALETA, CMEDIA_ALETA,
                              CRAIZ_ALETA, SW_ALETA)


# ----------------CARACTERÍSTICAS GEOMÉTRICAS DEL MISIL----------------
# Ángulo del cono (deg).
ANGULO_CONO = degrees(atan(.5 * DIAMETRO_M / LONGITUD_CONO))

# Superficie exterior del cono (m2)
SUP_CONO = pi * DIAMETRO_M / 2 * sqrt(LONGITUD_CONO**2 + DIAMETRO_M**2 / 4)

# Superficie de referencia del misil (m2)
SREF_MISIL = pi * DIAMETRO_M**2 / 4

SUP_TOTAL = pi * DIAMETRO_M * (LONGITUD_MISIL - LONGITUD_CONO)
# Superficie exterior del cilindro (m2)
SGASES = pi * (DIAMETRO_M * .45)**2
# Área de salida de los gases (consideramos el área de salida de la
# tobera, m2)
RATIO_AREAS = 1 - SGASES / SREF_MISIL  # Relación de áreas
TAO_ALETA = ESPESOR_ALETA / CMEDIA_ALETA  # TAO de la aleta
SWTOTAL_ALETAS = SW_ALETA * NUM_ALETAS  # Superficie total de aletas (m2)


def coef_resistencia_base_misil(mach_misil):
    '''Esta función define la resistencia de base del misil en función del
    Mach.
    '''

    if mach_misil < 0.8:
        return 0
    elif mach_misil < 1:
        term_indep = -1.548523
        term_uno = 6.05972764
        term_dos = -7.30548391
        term_tres = 2.96129532
        term_cuatro = 0
    elif mach_misil < 1.1:
        term_indep = 5.79090984e3
        term_uno = -2.19843314e4
        term_dos = 3.12774812e4
        term_tres = -1.97644892e4
        term_cuatro = 4.68059822e3
    elif mach_misil < 1.5:
        term_indep = -4.11856506
        term_uno = 1.42267421e1
        term_dos = -1.69678524e1
        term_tres = 8.771665
        term_cuatro = -1.67398037
    elif mach_misil < 2.2:
        term_indep = 3.0748e-1
        term_uno = -1.3258e-1
        term_dos = 2.8812e-2
        term_tres = 0
        term_cuatro = 0
    elif mach_misil <= 3.5:
        term_indep = 1.8481e-1
        term_uno = -2.2895e-2
        term_dos = 5.1876e-3
        term_tres = -4.0742e-4
        term_cuatro = 0
    elif mach_misil > 3.5:
        return .15

    return (term_cuatro * mach_misil**4 + term_tres * mach_misil**3 +
            term_dos * mach_misil**2 + term_uno * mach_misil +
            term_indep)


def cfcono_misil(re_cono, machl):
    '''Coeficiente de fricción del cono.
    '''
    # LAMINAR
    if re_cono < 1e6:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE
        cfi_cono = .664 * re_cono**(-1 / 2)
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE
        cf_cono = 2 * cfi_cono
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO
        cfm_cono = cf_cono / (1 + .17 * machl**2)**.1295
    # TURBULENTO
    else:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE
        cfi_cono = .288 / log10(re_cono)**2.45
        # CALCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE
        cf_cono = cfi_cono * 1.597 / log10(re_cono)**.15
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO
        cfm_cono = cf_cono / (1 + (GAMMA - 1) / 2 * machl**2)**.467
    return cfm_cono * SUP_CONO / SREF_MISIL


def cfcil(re_cilindro, machl):
    '''
    Coeficiente de fricción del cilindro.
    '''
    # LAMINAR
    if re_cilindro < 1e6:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE
        cfi_cil = .664 * re_cilindro**(-1 / 2)
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE
        cf_cil = 2 * cfi_cil
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO
        cfm_cil = cf_cil / (1 + .17 * machl**2)**.1295
    # TURBULENTO
    else:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE
        cfi_cil = .288 / log10(re_cilindro)**2.45
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE
        cf_cil = cfi_cil * 1.597 / log10(re_cilindro)**.15
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO
        cfm_cil = cf_cil / (1 + (GAMMA - 1) / 2 * machl**2)**.467
    return cfm_cil * SUP_TOTAL / SREF_MISIL


def cd_wave(mach, angulo, cd_f):
    '''Coeficiente de onda del misil.
    '''
    if mach >= 1:
        # RÉGIMEN SUPERSÓNICO.
        return (.083 + .096 / mach**2) * (angulo / 10)**1.69
    # RÉGIMEN SUBSÓNICO.
    ratio = LONGITUD_CONO / DIAMETRO_M
    return (60 / ratio**3 + .0025 * ratio) * cd_f


def cd_wave_aletas(mach):
    '''Coeficiente de onda de las aletas.
    '''
    # RÉGIMEN SUPERSÓNICO
    if mach >= 1:
        return 4 * TAO_ALETA**2 / (mach**2 - 1)**.5 * (SWTOTAL_ALETAS
                                                       / SREF_MISIL)
    # RÉGIMEN SUBSÓNICO.
    return 0


def cf_aletas(reyn_aleta, mach):
    '''Coeficiente de fricción de las aletas.
    '''
    # LAMINAR.
    if reyn_aleta < 1e6:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfialetas = .664 / reyn_aleta**.5
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL MEDIO.
        cf1aletas = 2 * cfialetas
        # CÁLCULO COEFICIENTE DE FRICCIÓN COMPRESIBLE.
        cfmaletas = cf1aletas / (1 + .17 * mach**2)**.1295
    # TURBULENTO.
    else:
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL INCOMPRESIBLE.
        cfialetas = .288 * (log10(reyn_aleta))**(-2.45)
        # CÁLCULO COEFICIENTE DE FRICCIÓN LOCAL COMPRESIBLE.
        cf1aletas = cfialetas * 1.597 * ((log10(reyn_aleta))**(-.15))
        # CÁLCULO COEFICIENTE DE FRICCIÓN MEDIO.
        cfmaletas = cf1aletas / (1 + (GAMMA - 1) / 2 * mach**2)**.467
    return cfmaletas * SWTOTAL_ALETAS / SREF_MISIL


def cdll(machl, alt):
    '''Coeficiente de resistencia total del misil.
    '''
    vel = machl * sqrt(GAMMA * R_AIR * temperature(alt))
    # CÁLCULO DEL COEFICIENTE DE RESISTENCIA BASE.
    cd_base_misil = coef_resistencia_base_misil(machl) * RATIO_AREAS
    re_cono = density(alt) * vel * LONGITUD_CONO / viscosity(alt)
    # Número de Reynolds en el cono.
    re_cil = density(alt) * vel * (LONGITUD_MISIL
                                   - LONGITUD_CONO) / viscosity(alt)
    # Número de Reynolds en el cilindro.
    # CÁLCULO DEL COEFICIENTE DE FRICCIÓN TOTAL REFERIDO A LA SUPERFICIE
    # TRANSVERSAL.
    cd_friccion_cono = cfcono_misil(re_cono, machl)
    # Coeficiente de fricción en el cono.
    cd_friccion_cil = cfcil(re_cil, machl)
    # Coeficiente de fricción en el cilindro.
    cd_friccion = cd_friccion_cono + cd_friccion_cil
    # CÁLCULO DEL COEFICIENTE DE ONDA.
    cd_onda = cd_wave(machl, ANGULO_CONO, cd_friccion)
    # RESISTENCIA DE LAS ALETAS.
    # COEFICIENTE DE ONDA.
    cd_onda_aletas = cd_wave_aletas(machl)
    re_aletas = density(alt) * vel * CRAIZ_ALETA / viscosity(alt)
    # Número de Reynolds en las aletas.
    # COEFICIENTE DE FRICCIÓN DE LAS ALETAS.
    cdfriccion_aletas = cf_aletas(re_aletas, machl)
    return (cd_base_misil + cd_friccion + cd_onda + cd_onda_aletas
            + cdfriccion_aletas)
