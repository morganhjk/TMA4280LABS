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

double ret (int n)
{
	return (sqrt (integral (n, &zeta) * 6.0));
}

void utest ()
{
	int n = 3;
	double expected = 2.85773803324704145;
	double computed = 0.0;
	char *message;

	computed = ret (n);

	if (expected == computed)
		message = "OK";
	else
		message = "FAIL";

	printf ("zeta0 utest: expected=%.17f, computed=%.17f, test=%s\n", expected, computed, message);
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
