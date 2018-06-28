# include <math.h>
# include "gravedad.h"
# include "modelo_msise00.h"
# include "velocidad_rotacional1.h"
# include "aero_misil.h"

double G0 = 9.81;

double *empuje(const double& gasto, const double& impulso, const double(&vel)[3])
{
	double empu[3];
	for (int i = 0; i < 3; ++i) {
		empu[i] = gasto * G0 * impulso * vel[i] / norm(vel);
	};
	double * pointer_empuje;
	pointer_empuje = &empu[0];
	return pointer_empuje;
}

double numero_mach(const double(&pos)[3], const double(&vel)[3])
{
	double altur = norm(pos) - RT;
	double tem = temperature(altur);
	double vel_sonido = sqrt(GAMMA * R_AIR * tem);
	double vel_aire[3], vel_relativa[3];
	double * c = cross(OMEGA_R, pos);
	for (int i = 0; i < 3; ++i) {
		vel_aire[i] = *(c + i);
		vel_relativa[i] = vel_aire[i] - vel[i];
	};
	return norm(vel_relativa) / vel_sonido;
	
}

double *resistencia(const double(&pos)[3], const double(&vel)[3])
{
	double altur = norm(pos) - RT;
	double mach = numero_mach(pos, vel);
	double cd_misil = cdll(mach, altur);
	double fuerza = 0.5 * GAMMA * pressure(altur) * pow(mach, 2) * SREF_MISIL * cd_misil;
	double resi[3];
	for (int i = 0; i < 3; ++i) {
		resi[i] = -fuerza * vel[i] / norm(vel);
	};
	double * pointer_resistencia;
	pointer_resistencia = &resi[0];
	return pointer_resistencia;
}

double *peso(const double(&pos)[3], const double& mas)
{
	double altur = norm(pos) - RT;
	double fuerza = mas * gravity(altur);
	double pe[3];
	for (int i = 0; i < 3; ++i) {
		pe[i] = -fuerza * pos[i] / (altur + RT);
	};
	double * pointer_peso;
	pointer_peso = &pe[0];
	return pointer_peso;
}

double aceleracion(int index, const double(&pos)[3], const double(&vel)[3], const double& mas, const double& gasto, const double& isp)
{
	double emp[3], res[3], pes[3], fuerza[3], acc[3];
	double * e = empuje(gasto, isp, vel);
	emp[0] = *e; emp[1] = *(e + 1); emp[2] = *(e + 2);
	double * r = resistencia(pos, vel);
	res[0] = *r; res[1] = *(r + 1); res[2] = *(r + 2);
	double * p = peso(pos, mas);
	pes[0] = *p; pes[1] = *(p + 1); pes[2] = *(p + 2);
	for (int i = 0; i < 3; ++i) {
		fuerza[i] = emp[i] + res[i] + pes[i];
		acc[i] = fuerza[i] / mas;
	};
	return acc[0];
}