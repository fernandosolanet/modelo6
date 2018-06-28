# include <iostream>
# include <string>
# include <stdio.h>
# include "integracion.h"
using namespace std;

int main()
{
	double t = 0.05;
	double x[3] = { 6390136.30, 20.5359212, 0.0 };
	double v[3] = { 409.49159427, 433.72234943, 0.0 };
	string m = "No imprimir";
	double dt = 0.1;
	double alt = 0;
	double mach = 1.5;
	double reynolds = 6e7;
	double mass = 20;
	double flowrate = 20;
	double specific_impulse = 300;
	double fuel_mass = 100;
	double total_mass = fuel_mass + mass;
	double maximum_altitude = 7.9e5;
	bool loss = false;

	printf("Masa: %.1f kg\n", mass);
	printf("Tiempo: %.3f s\n", t);
	printf("Posici\242n: (%.1f, %.1f, %.1f) m\n", x[0], x[1], x[2]);
	printf("Velocidad: (%.1f, %.1f, %.1f) m/s\n", v[0], v[1], v[2]);
	
	step(total_mass, t, x, v, flowrate, specific_impulse);
	// etapa(fuel_mass, total_mass, flowrate, specific_impulse, x, v, t, 0, dt, maximum_altitude, loss, m);

	printf("\n\nMasa: %.1f kg\n", total_mass);
	printf("Tiempo: %.3f s\n", t);
	printf("Posici\242n: (%.6g, %.6g, %.6g) m\n", x[0], x[1], x[2]);
	printf("Velocidad: (%.6g, %.6g, %.6g) m/s\n", v[0], v[1], v[2]);

	cout << "\n\nPress <enter> to exit";
	getline(cin, m);
}
