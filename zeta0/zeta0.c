#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
	S = pi^2 / 6
	S * 6 = pi^2
	sqrt (6S) = pi
*/

double zeta (int i)
{
	double di = (double) i;
	return 1.0 / (di * di);
}

double integral (int n, double (*f)(int))
{
	double accum = 0.0;
	
	for (int i = 1; i <= n; i++)
		accum += f(i);
	
	return accum;
}

int main (int argc, char **argv)
{
	if (argc != 2)
	{
		printf ("usage: ./zeta0 n\n");
		return 1;
	}
	
	int n = atoi (argv[1]);
	
	if (n <= 0)
	{
		printf ("n is bullshit, try again with different n\n");
		return 2;
	}
	
	printf ("Using n = %i\n", n);
	
	double ret = integral (n, &zeta);
	
	printf ("%.17f\n%.17f\n", ret, sqrt (ret * 6.0));
	
	return 0;
}
