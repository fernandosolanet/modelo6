#define _USE_MATH_DEFINES
#include <math.h>
#include "inputs_iniciales.h"

double SUP_TOTAL = M_PI * DIAMETRO_M * (LONGITUD_MISIL - LONGITUD_CONO);
double SUP_CONO = M_PI_2 * DIAMETRO_M * sqrt(pow(LONGITUD_CONO, 2) + pow(DIAMETRO_M / 2, 2));
double SREF_MISIL = M_PI_4*pow(DIAMETRO_M, 2);
double TAO_ALETA = ESPESOR_ALETA / CMEDIA_ALETA;
double SWTOTAL_ALETAS = SW_ALETA * NUM_ALETAS;
double SGASES = M_PI * pow(DIAMETRO_M * 0.45, 2);
double RATIO_AREAS = 1 - SGASES / SREF_MISIL;
double ANGULO_CONO = atan(0.5 * DIAMETRO_M / LONGITUD_CONO) * 180 / M_PI;

double coef_resistencia_base_misil(const double& mach_misil)
{
	double term_indep = 0;
	double term_uno = 0;
	double term_dos = 0;
	double term_tres = 0;
	double term_cuatro = 0;
	if (mach_misil < 0.8) {
		return 0;
	}
	else if (mach_misil < 1) {
		term_indep = -1.548523;
		term_uno = 6.05972764;
		term_dos = -7.30548391;
		term_tres = 2.96129532;
	}
	else if (mach_misil < 1.1) {
		term_indep = 5.79090984e3;
		term_uno = -2.19843314e4;
		term_dos = 3.12774812e4;
		term_tres = -1.97644892e4;
		term_cuatro = 4.68059822e3;
	}
	else if (mach_misil < 1.5) {
		term_indep = -4.11856506;
		term_uno = 1.42267421e1;
		term_dos = -1.69678524e1;
		term_tres = 8.771665;
		term_cuatro = -1.67398037;
	}
	else if (mach_misil < 2.2) {
		term_indep = 3.0748e-1;
		term_uno = -1.3258e-1;
		term_dos = 2.8812e-2;
	}
	else if (mach_misil <= 3.5) {
		term_indep = 1.8481e-1;
		term_uno = -2.2895e-2;
		term_dos = 5.1876e-3;
		term_tres = -4.0742e-4;
	}
	else if (mach_misil > 3.5) {
		term_indep = 0.15;
	};
	return term_cuatro * pow(mach_misil, 4) + term_tres * pow(mach_misil, 3) + term_dos * pow(mach_misil, 2) + term_uno * mach_misil + term_indep;
}

double cfcono_misil(const double& re_cono, const double& machl)
{
	double cfi_cono, cf_cono, cfm_cono;
	if (re_cono < 1e6) {
		cfi_cono = 0.664 / sqrt(re_cono);
		cf_cono = 2 * cfi_cono;
		cfm_cono = cf_cono * pow((1 + 0.17 * pow(machl, 2)), -0.1295);
	}
	else {
		cfi_cono = 0.288 * pow(log10(re_cono), -2.45);
		cf_cono = cfi_cono * 1.597 * pow(log10(re_cono), -0.15);
		cfm_cono = cf_cono * pow((1 + (GAMMA - 1) / 2 * pow(machl, 2)), -0.467);
	};
	return cfm_cono * SUP_CONO / SREF_MISIL;
}

double cfcil(const double& re_cilindro, const double& machl)
{
	double cfi_cil, cf_cil, cfm_cil;
	if (re_cilindro < 1e6) {
		cfi_cil = 0.664 / sqrt(re_cilindro);
		cf_cil = 2 * cfi_cil;
		cfm_cil = cf_cil * pow((1 + 0.17 * pow(machl, 2)), -0.1295);
	}
	else {
		cfi_cil = 0.288 * pow(log10(re_cilindro), -2.45);
		cf_cil = cfi_cil * 1.597 * pow(log10(re_cilindro), -0.15);
		cfm_cil = cf_cil * pow((1 + (GAMMA - 1) / 2 * pow(machl, 2)), -0.467);
	};
	return cfm_cil * SUP_TOTAL / SREF_MISIL;
}

double cd_wave(const double& mach, const double& angulo, const double& cd_f)
{
	if (mach >= 1) {
		return (0.083 + 0.096 * pow(mach, -2)) * pow(angulo / 10, 1.69);
	};
	double ratio = LONGITUD_CONO / DIAMETRO_M;
	return (60 / pow(ratio, 3) + 0.0025 * ratio) * cd_f;
}

double cd_wave_aletas(const double& mach)
{
	if (mach < 1) {
		return 0;
	}
	else {
		return 4 * pow(TAO_ALETA, 2) / sqrt(pow(mach, 2) - 1) * SWTOTAL_ALETAS / SREF_MISIL;
	};
}


double cf_aletas(const double& reyn_aleta, const double& mach)
{
	double cfialetas, cf1aletas, cfmaletas;
	if (reyn_aleta < 1e6) {
		cfialetas = 0.664 / sqrt(reyn_aleta);
		cf1aletas = 2 * cfialetas;
		cfmaletas = cf1aletas * pow(1 + 0.17 * pow(mach, 2), -0.1295);
	}	
	else {
		cfialetas = 0.288 * pow(log10(reyn_aleta), -2.45);
		cf1aletas = cfialetas * 1.597 * pow(log10(reyn_aleta), -0.15);
		cfmaletas = cf1aletas * pow(1 + (GAMMA - 1) / 2 * pow(mach, 2), -0.467);
	};
	return cfmaletas * SWTOTAL_ALETAS / SREF_MISIL;
}

double cdll(const double& machl, const double& alt)
{
	double vel = machl * sqrt(GAMMA * R_AIR * temperature(alt));
	double cd_base_misil = coef_resistencia_base_misil(machl) * RATIO_AREAS;
	double re_cono = density(alt) * vel * LONGITUD_CONO / viscosity(alt);
	double re_cil = density(alt) * vel * (LONGITUD_MISIL - LONGITUD_CONO) / viscosity(alt);
	double cd_friccion_cono = cfcono_misil(re_cono, machl);
	double cd_friccion_cil = cfcil(re_cil, machl);
	double cd_friccion = cd_friccion_cono + cd_friccion_cil;
	double cd_onda = cd_wave(machl, ANGULO_CONO, cd_friccion);
	double cd_onda_aletas = cd_wave_aletas(machl);
	double re_aletas = density(alt) * vel * CRAIZ_ALETA / viscosity(alt);
	double cdfriccion_aletas = cf_aletas(re_aletas, machl);
	return cd_base_misil + cd_friccion + cd_onda + cd_onda_aletas + cdfriccion_aletas;
}


