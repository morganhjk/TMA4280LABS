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
	return sqrt (integral (n, &zeta) * 6.0);
}

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int utest (void)
{
	// Test name
	char *test_name = "zeta0 utest";

	// Set expected value and get computed value
	double expected = 2.85773803324704145368;
	double computed = ret (3);

	// Set test message
	char *message = (expected == computed) ? "OK" : "FAIL";

	// Print results
	printf ("%s: expected=%.20f, computed=%.20f, test=%s\n",
		test_name, expected, computed, message);

	// Done
	return 0;
}

int vtest (void)
{
	// Test name and paths to log file
	char *test_name = "zeta0 vtest";
	char *log_rel_path = "vtest.txt";
	char *log_abs_path;

	// Open log file for writing
	FILE *log = fopen (log_rel_path, "w");

	// Main loop
	for (int i = 1; i <= 24; i++)
	{
		// Set n to a power of 2
		int n = 2 << i;

		// Get computed value and calculate error
		double computed = ret (n);
		double error = fabs (M_PI - computed);

		// Write to log
		fprintf (log, "%s: computed=%.20f, error=%.20f, n=%i\n",
			test_name, computed, error, n);
	}

	// Close log file
	fclose (log);

	// Get absolute path to log file and print
	log_abs_path = realpath (log_rel_path, NULL);
	printf ("%s results written to: %s\n", test_name, log_abs_path);
	free (log_abs_path);

	// Done
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

int main (__attribute__ ((unused)) int argc, __attribute__ ((unused)) char **argv)
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

	// Get number of iterations (n)
	int n = atoi (argv[1]);

	// Check n
	if (n <= 0)
	{
		fprintf (stderr, "'%i' is not a positive integer, try again with a different <n>\n", n);
		usage (argv[0]);
		return 1;
	}

	// Print and exit
	printf ("%.20f\n", ret (n));
	return 0;
#endif
}
