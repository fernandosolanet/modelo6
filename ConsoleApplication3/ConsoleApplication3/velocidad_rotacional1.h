# define _USE_MATH_DEFINES
# include <math.h>
double pi = 3.141592653589793238462643383279502884;

double ANO_MEDIO = (265 * 400 + 100 - 1) / 400;
double PERIODO = ANO_MEDIO / (ANO_MEDIO + 1) * 24 * 60 * 60;
double V_ANGULAR = 2 * M_PI / PERIODO;
double OMEGA_R[3] = { 0, 0, V_ANGULAR };
