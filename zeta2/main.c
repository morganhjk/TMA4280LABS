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

// Vector
double *numbers;

///////////////////////////////////////////////////////////////////////////////
// Math functions
///////////////////////////////////////////////////////////////////////////////
double zeta (int i)
{
	double di = (double) i;
	return 1.0 / (di * di);
}

double integral (int from, int to, double (*f)(int))
{
	double accum = 0.0;

	if (from && to)
#pragma omp parallel for num_threads(nthreads) schedule(static) reduction(+:accum)
		for (int i = from; i <= to; i++)
			accum += f(i);

	return accum;
}

double ret (void)
{
	double accum = 0.0;

	for (int i = 0; i < commsize; i++)
		accum += numbers[i];

	return sqrt (accum * 6.0);
}

///////////////////////////////////////////////////////////////////////////////
// MPI functions
///////////////////////////////////////////////////////////////////////////////
void init (int n)
{
	// Allocate global indexes
	gi = malloc (commsize * 2 * sizeof (int));

	// Allocate vector of elements to sum
	numbers = malloc (commsize * sizeof (double));

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

void worker (void)
{
	// Scatter 'from' and 'to' indexes to all ranks
	MPI_Scatter (gi, 2, MPI_INT, li, 2, MPI_INT, 0, MPI_COMM_WORLD);

	// Calculate my range of elements (local sum)
	double lsum = integral (li[0], li[1], &zeta);

	// Gather sums from all ranks
	MPI_Gather (&lsum, 1, MPI_DOUBLE, numbers, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void clean (void)
{
	// Free globals
	free (numbers);
	free (gi);
}

double avgtime (double start, double end)
{
	// Declare array for all walltimes
	double times[commsize];

	// Calculate my walltime
	double time = end - start;

	// Gather walltimes from all ranks
	MPI_Gather (&time, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Calculate average walltime
	double avg = 0.0;
	if (!myrank)
	{
		for (int i = 0; i < commsize; i++)
			avg += times[i];

		avg /= commsize;
	}

	return avg;
}

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int vtest (void)
{
	// Test name and paths to log file
	char *test_name = "zeta2 vtest";
	char *log_rel_path = "vtest.txt";
	char *log_abs_path;

	// Pointer to log file
	FILE *log = NULL;

	// Helper table for LaTeX
#ifdef TEXTABLE
	double errtable[3][24];
	double timtable[3][24];
#endif

	// Open log file for writing (append mode)
	if (!myrank)
		log = fopen (log_rel_path, "a");

	// Main loop
	for (int i = 1; i <= 3; i++)
	{
		nthreads = pow (2, i);

		for (int j = 1; j <= 24; j++)
		{
			// Start time
			double t1 = MPI_Wtime ();

			// Set n to a power of 2
			int n = 2 << j;

			// Allocate globals and generate indexes
			if (!myrank)
				init (n);

			// Perform work
			worker ();

			// End time
			double t2 = MPI_Wtime ();

			// Average time across procs (returns sensible result only on rank 0)
			double tavg = avgtime (t1, t2);

			// Get computed value and calculate error, then write to log
			if (!myrank)
			{
				double computed = ret ();
				double error = fabs (M_PI - computed);

#ifdef TEXTABLE
            	errtable[i-1][j-1] = error;
            	timtable[i-1][j-1] = tavg;
#else
				fprintf (log, "%s p=%i t=%i: computed=%.20f, error=%.20f, time=%.20f, n=%i\n",
						 test_name, commsize, nthreads, computed, error, tavg, n);
#endif
			}
		}
	}

#ifdef TEXTABLE
	if (!myrank)
	{
		for (int j = 1; j <= 24; j++)
			fprintf (log, "\t\t%i\t& %.20f & %.20f & %.20f \\\\\n",
				2 << j, errtable[0][j-1], errtable[1][j-1], errtable[2][j-1]);

		for (int j = 1; j <= 24; j++)
			fprintf (log, "\t\t%i\t& %.20f & %.20f & %.20f \\\\\n",
				2 << j, timtable[0][j-1], timtable[1][j-1], timtable[2][j-1]);
	}
#endif

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

	// Get number of iterations (n)
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

	// Get number of threads (nthreads)
	nthreads = atoi (argv[2]);

	// Check nthreads
	if (nthreads <= 0)
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a positive integer, try again with a different <t>\n", nthreads);
			usage (argv[0]);
		}
		return 3;
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
