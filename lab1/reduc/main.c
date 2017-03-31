///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////
// Global variables
///////////////////////////////////////////////////////////////////////////////
// Size and rank
int commsize;
int myrank;

// From/to indexes
int *gi;	// Global
int li[2];	// Local

// Vector
double *numbers;

// Local calculation
double mycalculation;

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
		for (int i = from; i <= to; i++)
			accum += f(i);

	return accum;
}

double ret (double sum)
{
	return sqrt (sum * 6.0);
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
	mycalculation = integral (li[0], li[1], &zeta);
}

double reducesum (void)
{
#ifdef MPISUM
	// Calculate sum using allreduce
	double sum = 0.0;
	MPI_Allreduce (&mycalculation, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sum;
#elif RECURSIVE
	// Variables
	double sum = mycalculation;
	double recdat = 0.0;
	int sendto;
	MPI_Status stat;

	// Run for log n steps
	for (int step = 1; step < commsize; step = step << 1)
	{
		// If we're odd, exchange with previous, if not, exchange with next
		if ((myrank / step) & 0x01)
			sendto = myrank - step;
		else
			sendto = myrank + step;

		// Wrap-around
		sendto = sendto % commsize;

		// Exchange
		MPI_Sendrecv (&sum, 1, MPI_DOUBLE, sendto, 0,
						&recdat, 1, MPI_DOUBLE, sendto, 0,
						MPI_COMM_WORLD, &stat);

		// Add to sum
		sum += recdat;
	}

	// Done
	return sum;
#else
	// Gather sums from all ranks
	MPI_Gather (&mycalculation, 1, MPI_DOUBLE, numbers, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Sum vector elements
	double sum = 0.0;

	if (!myrank)
	{
		for (int i = 0; i < commsize; i++)
			sum += numbers[i];
	}

	return sum;
#endif
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
// Main functions
///////////////////////////////////////////////////////////////////////////////
void usage (char *appname)
{
	printf ("usage: mpirun -np <p> %s <n>\n", appname);
	printf ("where:\n");
	printf (" <p> is a power of two (number of processes)\n");
	printf (" <n> is a positive integer (number of iterations)\n");
}

int main (int argc, char **argv)
{
	// MPI Init
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Check number of arguments
	if (argc != 2)
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

	// Start time
	double t1 = MPI_Wtime ();

	// Allocate globals and generate indexes
	if (!myrank)
		init (n);

	// Perform work
	worker ();

	// Reduce each process partial sum to one
	double sum = reducesum ();

	// End time
	double t2 = MPI_Wtime ();

	// Average time across procs (returns sensible result only on rank 0)
	double tavg = avgtime (t1, t2);

	// Final work, print, and cleanup
	if (!myrank)
	{
		printf ("computed=%.20f, time=%.20f\n", ret (sum), tavg);
		clean ();
	}

	// Done
	MPI_Finalize ();
	return 0;
}
