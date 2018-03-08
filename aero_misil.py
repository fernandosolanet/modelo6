# -*- coding: utf-8 -*-
"""

@author: Team REOS

Funciones que nos dan los coeficientes aerodinámicos del avión F4.
Coeficientes aerodinámicos de sustentación y de resistencia y todo aquello
relacionado con los mismos.

"""

from math import log10, pi, degrees, atan
from modeloISA import GAMMA, density, temperature, R_AIR, viscosity


#--------------------CARACTERÍSTICAS GEOMÉTRICAS DEL MISIL--------------------

DIAMETRO_M = .5  # Diámetro del misil (m).
LONGITUD_CONO = .9  # Longitud del cono del misil (m).
LONGITUD_MISIL = 3  # Longitud total del misil (m).
ANGULO_CONO = degrees(atan(.5 * DIAMETRO_M / LONGITUD_CONO))
# Ángulo del cono (deg).
SUP_CONO = pi * DIAMETRO_M / 2 * (LONGITUD_CONO**2 + DIAMETRO_M**2 / 4)**.5
# Superficie exterior del cono (m2).
SREF_MISIL = pi * DIAMETRO_M**2 / 4  # Superficie de referencia del misil (m2).
SUP_TOTAL = pi * DIAMETRO_M * (LONGITUD_MISIL - LONGITUD_CONO)
# Superficie exterior del misil (m2).
SGASES = pi * (DIAMETRO_M * .45)**2
# Área de salida de los gases (consideramos el área de salida de la tobera,
# m2).
RATIO_AREAS = 1 - SGASES / SREF_MISIL  # Relación de áreas.

ESPESOR_ALETA = .0065  # Espesor de la aleta (m).
CMEDIA_ALETA = .18  # Cuerda media de la aleta (m).
CRAIZ_ALETA = .24  # Cuerda raiz de la aleta (m).
TAO_ALETA = ESPESOR_ALETA / CMEDIA_ALETA  # TAO de la aleta.
SW_ALETA = .07875  # Superficie de una aleta del AIM (tomado como ref., m2).
NUM_ALETAS = 4  # Número de aletas.
SWTOTAL_ALETAS = SW_ALETA * NUM_ALETAS  # Superficie total de aletas (m2).

def coef_resistencia_base_misil(Ml):

	if Ml < 0.8:
		x0 = 0
		x1 = 0
		x2 = 0
		x3 = 0
		x4 = 0
	elif Ml < 1:
		x0 = -1.548523
		x1 = 6.05972764
		x2 = -7.30548391
		x3 = 2.96129532
		x4 = 0
	elif Ml < 1.1:
		x0 = 5.79090984*10**3
		x1 = -2.19843314*10**4
		x2 = 3.12774812*10**4
		x3 = -1.97644892*10**4
		x4 = 4.68059822*10**3
	elif Ml < 1.5:
		x0 = -4.11856506
		x1 = 1.42267421*10**1
		x2 = -1.69678524*10**1
		x3 = 8.771665
		x4 = -1.67398037
	elif Ml < 2.2:
		x0 = 3.0748*10**-1
		x1 = -1.3258*10**-1
		x2 = 2.8812*10**-2
		x3 = 0
		x4 = 0
	elif Ml <=3.5:                            
		x0 = 1.8481*10**-1
		x1 = -2.2895*10**-2
		x2 = 5.1876*10**-3
		x3 = -4.0742*10**-4
		x4 = 0
	elif Ml >3.5:                            
		x0 = 0.15
		x1 = 0
		x2 = 0
		x3 = 0
		x4 = 0
		    		
	return x4*Ml**4 + x3*Ml**3 + x2*Ml**2 + x1*Ml + x0

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

def cdll(machl, alt):
    '''Coeficiente de resistencia total del misil.
    '''
    vel = machl * (GAMMA * R_AIR * temperature(alt))**.5
    #CÁLCULO DEL COEFICIENTE DE RESISTENCIA BASE.
    cd_base_misil = coef_resistencia_base_misil(machl) * RATIO_AREAS
    re_cono = density(alt) * vel * LONGITUD_CONO / viscosity(alt)
    # Número de Reynolds en el cono.
    re_cil = density(alt) * vel * (LONGITUD_MISIL
                                   - LONGITUD_CONO) / viscosity(alt)
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
    re_aletas = density(alt) * vel * CRAIZ_ALETA / viscosity(alt)
    # Número de Reynolds en las aletas.
    #COEFICIENTE DE FRICCIÓN DE LAS ALETAS.
    cdfriccion_aletas = cf_aletas(re_aletas, machl)
    return (cd_base_misil + cd_friccion + cd_onda + cd_onda_aletas
            + cdfriccion_aletas)
