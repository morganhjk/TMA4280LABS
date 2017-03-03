///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

///////////////////////////////////////////////////////////////////////////////
// Global variables
///////////////////////////////////////////////////////////////////////////////
// Size and rank
int commsize;
int myrank;

// Number of threads per process
int nthreads = 1;

// From/to indexes
int *gi;	// Global
int li[2];	// Local

// Vectors
double *numbersa;
double *numbersb;

///////////////////////////////////////////////////////////////////////////////
// Math functions
///////////////////////////////////////////////////////////////////////////////
double machin (int i, double x)
{
	double ai = (double) ((2 * i) - 1);
	double frac = pow (x, ai) / ai;
	return ((i-1) % 2) ? (-1.0 * frac) : frac;
}

double integral (int from, int to, double x, double (*f)(int, double))
{
	double accum = 0.0;

	if (from && to)
#pragma omp parallel for num_threads(nthreads) schedule(static) reduction(+:accum)
		for (int i = from; i <= to; i++)
			accum += f(i, x);

	return accum;
}

double ret (void)
{
	double accuma = 0.0;
	double accumb = 0.0;

	for (int i = 0; i < commsize; i++)
	{
		accuma += numbersa[i];
		accumb += numbersb[i];
	}

	return 4.0 * (4.0 * accuma - accumb);
}

///////////////////////////////////////////////////////////////////////////////
// MPI functions
///////////////////////////////////////////////////////////////////////////////
void init (int n)
{
	// Allocate global indexes
	gi = malloc (commsize * 2 * sizeof (int));

	// Allocate vectors of elements to sum
	numbersa = malloc (commsize * sizeof (double));
	numbersb = malloc (commsize * sizeof (double));

	// Divide the number of iterations on the number of procs
	// to get the number of equal-sized parts
	int parts = n / commsize;

	for (int r = 0, s = 0; r < commsize; r++, s+=2)
	{
		// Set 'from' index for the current rank
		gi[s] = parts * r + 1;

		// Set 'to' index for the current rank
		if (r != commsize-1)
			gi[s+1] = parts * (r + 1);
		else
		{
			// Last rank gets the remaining iterations
			gi[s+1] = n;
		}
	}
}

void worker ()
{
	// Scatter 'from' and 'to' indexes to all ranks
	MPI_Scatter (gi, 2, MPI_INT, li, 2, MPI_INT, 0, MPI_COMM_WORLD);

	// Calculate my range of elements (local sums)
	double lsuma = integral (li[0], li[1], (1.0 / 5.0), &machin);
	double lsumb = integral (li[0], li[1], (1.0 / 239.0), &machin);

	// Gather sums from all ranks
	MPI_Gather (&lsuma, 1, MPI_DOUBLE, numbersa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather (&lsumb, 1, MPI_DOUBLE, numbersb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void clean (void)
{
	// Free globals
	free (numbersa);
	free (numbersb);
	free (gi);
}

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int vtest (void)
{
	// Test name and paths to log file
	char *test_name = "mach1 vtest";
	char *log_rel_path = "vtest.txt";
	char *log_abs_path;

	// Pointer to log file
	FILE *log = NULL;

	// Open log file for writing (append mode)
	if (!myrank)
		log = fopen (log_rel_path, "a");

	// Main loop
	for (int i = 1; i <= 24; i++)
	{
		// Set n to a power of 2
		int n = 2 << i;

		// Allocate globals and generate indexes
		if (!myrank)
			init (n);

		// Perform work
		worker ();

		// Get computed value and calculate error, then write to log
		if (!myrank)
		{
			double computed = ret ();
			double error = fabs (M_PI - computed);

			fprintf (log, "%s p=%i: computed=%.20f, error=%.20f, n=%i\n",
					 test_name, commsize, computed, error, n);
		}
	}

	if (!myrank)
	{
		// Close log file
		fclose (log);

		// Get absolute path to log file and print
		log_abs_path = realpath (log_rel_path, NULL);
		printf ("%s p=%i results written to: %s\n", test_name, commsize, log_abs_path);
		free (log_abs_path);
	}

	// Done
	MPI_Finalize ();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Main functions
///////////////////////////////////////////////////////////////////////////////
void usage (char *appname)
{
	printf ("usage: mpirun -np <p> %s <n> <t>\n", appname);
	printf ("where:\n");
	printf (" <p> is a power of two (number of processes)\n");
	printf (" <n> is a positive integer (number of iterations)\n");
	printf (" <t> is a positive integer (number of threads per process)\n");
}

int main (int argc, char **argv)
{
	// MPI Init
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

#ifdef VTEST
	// Check commsize
	if (((commsize & (commsize - 1)) != 0))
		return 1;

	// Perform verification test
	return vtest ();
#else
	// Check number of arguments
	if (argc != 3)
	{
		if (!myrank)
			usage (argv[0]);
		return 0;
	}

	// Check commsize
	if (((commsize & (commsize - 1)) != 0))
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a power of two, try again with a different <p>\n", commsize);
			usage (argv[0]);
		}
		return 1;
	}

	// Get n from arguments
	int n = atoi (argv[1]);

	// Check n
	if (n <= 0)
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a positive integer, try again with a different <n>\n", n);
			usage (argv[0]);
		}
		return 2;
	}

	// Get number of threads
	nthreads = atoi (argv[2]);

	// Check n
	if (nthreads <= 0)
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a positive integer, try again with a different <n>\n", nthreads);
			usage (argv[0]);
		}
		return 2;
	}

	// Allocate globals and generate indexes
	if (!myrank)
		init (n);

	// Perform work
	worker ();

	// Final work, print, and cleanup
	if (!myrank)
	{
		printf ("%.20f\n", ret ());
		clean ();
	}

	// Done
	MPI_Finalize ();
	return 0;
#endif
}
