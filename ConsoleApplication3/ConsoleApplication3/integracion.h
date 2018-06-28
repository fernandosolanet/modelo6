# include <limits>
# include <string>
# include <fstream>
# include <iostream>
# include "numpy.h"
# include "mecanica.h"

double DT = 0.05;

double step(double& mas, double& tie, double(&pos)[3], double(&vel)[3], const double& gasto, const double& isp, double vloss = 0.0, const double& masa_minima = 0, double& step_size = DT, const bool& perdidas = false)
{
	double dtl = step_size;
	double acel[3];
	double masa = mas - gasto * dtl;
	if (masa < masa_minima) {
		dtl = (mas - masa_minima) / gasto;
		masa = masa_minima;
	};
	tie = tie + dtl;
	double * acceleration;
	acel[0] = aceleracion(pos, vel, (mas + masa) / 2, gasto, isp);
	acceleration = &acel[0];
	printf("\n%g", *acceleration);
	printf("\n%g", *acceleration);
	printf("\n%g", *acceleration);
	acel[1] = *(acceleration + 1); acel[2] = *(acceleration + 2);
	printf("\nAceleraci\242n : (%.6g, %.6g, %.6g) m/s2", acel[0], acel[1], acel[2]);
	for (int i = 0; i < 3; ++i) {
		pos[i] = pos[i] + vel[i] * dtl + 0.5 * acel[i] * pow(dtl, 2);
		vel[i] = vel[i] + acel[i] * dtl;
	};
	if (perdidas) {
		double av_pos[3], av_vel[3], resis[3], weight[3], av_mas, loss_aero, loss_grav;
		double *resi, *wei;
		av_mas = (mas + masa) / 2;
		for (int i = 0; i < 3; ++i) {
			av_pos[i] = pos[i] - 0.5 * vel[i] * dtl - 0.25 * acel[i] * pow(dtl, 2);
			av_vel[i] = vel[i] - 0.5 * acel[i] * dtl;
		};
		resi = resistencia(av_pos, av_vel);
		wei = peso(av_pos, av_mas);
		for (int i = 0; i < 3; ++i) {
			resis[i] = *(resi + i);
			weight[i] = *(wei + i);
		};
		double resistance = norm(resis);
		loss_aero = dtl * resistance / av_mas;
		loss_grav = -dtl * dot(weight, av_vel) / (norm(av_vel) * av_mas);
		vloss = vloss + loss_aero + loss_grav;
	};
	mas = masa;
	return vloss;
}

double etapa(double& masa_combustible, double& masa_total, const double& gasto, const double& isp, double(&pos)[3], double(&vel)[3], double& t, double vloss = 0.0, double& step_size = DT, const double& altura_maxima = numeric_limits<double>::infinity(), const bool& perdidas = false, const string& imprimir = "No_imprimir")
{
	double resto = masa_total - masa_combustible;
	bool consumido = (masa_combustible <= 0);
	double altur = norm(pos) - RT;
	ofstream archivo(imprimir, ios::trunc);
	char output[100];
	while ( (altur < altura_maxima) && (dot(pos, vel) >= 0) && (!consumido) ){
		vloss = step(masa_total, t, pos, vel, gasto, isp, vloss, resto, step_size, perdidas);
		altur = norm(pos) - RT;
		masa_combustible = masa_total - resto;
		consumido = (masa_combustible <= 0);
		if (archivo.is_open()) {
			int n = sprintf_s(output, 100, "\n%.3f\t%.0f\t%.0f\t%.1f", t, altur, norm(vel), masa_total);
			archivo << output;
		};
	}
	return vloss;
}

double vuelo_libre(const double& masa, double(&pos)[3], double(&vel)[3], double& t, const double& t_de_vuelo = numeric_limits<double>::infinity(), double vloss = 0, double& step_size = DT, const double& altura_maxima = numeric_limits<double>::infinity(), const bool& perdidas = false, const string& imprimir = "No_imprimir")
{
	return 0.0;
}