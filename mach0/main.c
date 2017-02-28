///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
// Math functions
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int utest (void)
{
	double expected = 3.14162102932503462;
	double computed = ret (3);
	char   *message = (expected == computed) ? "OK" : "FAIL";

	printf ("mach0 utest: expected=%.17f, computed=%.17f, test=%s\n",
		expected, computed, message);

	return 0;
}

int vtest (void)
{
	FILE *fp = fopen ("test-results.txt", "w");

	for (int i = 1; i <= 24; i++)
	{
		int n = 2 << i;
		double pi = 3.14159265358979323;
		double computed = ret (n);
		double error = pi - computed;

		if (error < 0.0)
			error = -error;

		fprintf (fp, "mach0 vtest: computed=%.17f, error=%.17f, n=%i\n",
			computed, error, n);
	}

	fclose (fp);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main functions
///////////////////////////////////////////////////////////////////////////////
void usage (char *appname)
{
	printf ("usage: %s <n>\n", appname);
	printf ("where:\n");
	printf (" <n> is a positive integer (number of iterations)\n");
}

int main (int argc, char **argv)
{
#ifdef UTEST
	return utest ();
#elif VTEST
	return vtest ();
#else
	// Check number of arguments
	if (argc != 2)
	{
		usage (argv[0]);
		return 0;
	}

	// Get n from arguments
	int n = atoi (argv[1]);

	// Check n
	if (n <= 0)
	{
		fprintf (stderr, "'%i' is not a positive integer, try again with a different <n>\n", n);
		usage (argv[0]);
		return 1;
	}

	// Print and exit
	printf ("%.17f\n", ret (n));
	return 0;
#endif
}
