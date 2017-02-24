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

double ret (int n)
{
	double reta = integral (n, (1.0 / 5.0), &machin);
	double retb = integral (n, (1.0 / 239.0), &machin);

	return (4.0 * (4.0 * reta - retb));
}

void utest ()
{
	int n = 3;
	double expected = 3.14162102932503462;
	double computed = 0.0;
	char *message;

	computed = ret (n);

	if (expected == computed)
		message = "OK";
	else
		message = "FAIL";

	printf("mach0 utest: expected=%.17f, computed=%.17f, test=%s\n", expected, computed, message);
}

int main (int argc, char **argv)
{
#ifdef UTEST
	utest();
	return 0;
#endif

	if (argc != 2)
	{
		printf ("usage: %s n\n", argv[0]);
		return 1;
	}
	
	int n = atoi (argv[1]);
	
	if (n <= 0)
	{
		printf ("n is bullshit, try again with different n\n");
		return 2;
	}
	
	printf ("%.17f\n", ret (n));
	return 0;
}
