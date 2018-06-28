#include <math.h>
#include <stdio.h>
#include <numeric>
#include <string>
using namespace std;

double *cross(const double (&a)[3], const double (&b)[3])
{
	double c[3];
	double * d;
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	d = &c[0];
	return d;
}

double dot(const double (&a)[3], const double (&b)[3])
{
	double c;
	c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return c;
}

double norm(const double (&a)[3])
{
	double b = 0;
	for (const double& c : a)
	{
		b = b + c*c;
	}
	return sqrt(b);
}
