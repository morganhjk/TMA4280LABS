#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double machin (int i, double x)
{
	double ai = (double) ((2 * i) - 1);
	double frac = pow (x, ai) / ai;
	return ((i-1) % 2) ? (-1.0 * frac) : frac;
}

double integral (int n, double x, double (*f)(int, double))
{
	double accum = 0.0;
	
	for (int i = 1; i <= n; i++)
		accum += f(i, x);
	
	return accum;
}

int main (int argc, char **argv)
{
	if (argc != 2)
	{
		printf ("usage: ./mach0 n\n");
		return 1;
	}
	
	int n = atoi (argv[1]);
	
	if (n <= 0)
	{
		printf ("n is bullshit, try again with different n\n");
		return 2;
	}
	
	double reta = integral (n, (1.0 / 5.0), &machin);
	double retb = integral (n, (1.0 / 239.0), &machin);
	double ret = 4.0 * (4 * reta - retb);
	
	printf ("%.17f\n", ret);
	return 0;
}
