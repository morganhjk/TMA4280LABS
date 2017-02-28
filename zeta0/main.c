///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
// Math functions
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int utest (void)
{
	double expected = 2.85773803324704145;
	double computed = ret (3);
	char   *message = (expected == computed) ? "OK" : "FAIL";

	printf ("zeta0 utest: expected=%.17f, computed=%.17f, test=%s\n",
		expected, computed, message);

	return 0;
}

int vtest (void)
{
	for (int i = 1; i <= 24; i++)
	{
		int n = 2 << i;
		double pi = 3.14159265358979323;
		double computed = ret (n);

		printf ("zeta0 vtest: computed=%.17f, error=%.17f, n=%i\n",
			computed, (pi - computed), n);
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main functions
///////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv)
{
#ifdef UTEST
	return utest();
#elif VTEST
	return vtest ();
#else
	if (argc != 2)
	{
		printf ("usage: %s n\n", argv[0]);
		return 1;
	}

	int n = atoi (argv[1]);

	if (n <= 0)
	{
		printf ("n is bullshit, try again with different n (got %i)\n", n);
		return 2;
	}

	printf ("%.17f\n", ret (n));
	return 0;
#endif
}
