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
int nthreads;

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
void lille (void)
{
	// Test name and paths to log file
	char *test_name = "lille";
	char *log_rel_path, *log_abs_path;
	asprintf (&log_rel_path, "%s-p%i-t%i.txt", test_name, commsize, nthreads);

	// Pointer to log file
	FILE *log = NULL;

	// Open log file for writing
	if (!myrank)
		log = fopen (log_rel_path, "w");

	// Main loop
	for (int i = 1; i <= 25; i++)
	{
		// Start time
		double t1 = MPI_Wtime ();

		// Set n to a power of 2
		int n = 2 << i;

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

			fprintf (log, "%s p=%i t=%i: computed=%.20f, error=%.20f, time=%.20f, n=%i\n",
					 test_name, commsize, nthreads, computed, error, tavg, n);
		}
	}

	if (!myrank)
	{
		// Close log file
		fclose (log);

		// Get absolute path to log file and print
		log_abs_path = realpath (log_rel_path, NULL);
		printf ("%s (p=%i t=%i) results written to: %s\n", test_name, commsize, nthreads, log_abs_path);

		// Cleanup
		free (log_abs_path);
		free (log_rel_path);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Main functions
///////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv)
{
	// MPI Init
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// OpenMP Init
	nthreads = omp_get_num_threads ();

	// Ensure that number of procs is a power of 2
	if (((commsize & (commsize - 1)) != 0))
		return 1;

	// Perform test
	lille ();

	// Done
	MPI_Finalize ();
	return 0;
}
