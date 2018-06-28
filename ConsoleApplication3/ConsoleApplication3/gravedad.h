#include <math.h>
using namespace std;

double MT = 5.972e24;
double G = 6.673e-11;
double MU = G * MT;
double RT = 6378136.3;

double gravity(const double& alt)
{
	double radio = RT + alt;
	return MU / pow(radio, 2);
}
