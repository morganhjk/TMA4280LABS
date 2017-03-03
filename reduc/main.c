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

double integral (int from, int to, double (*f)(int))
{
	double accum = 0.0;
	
	if (from && to)
		for (int i = from; i <= to; i++)
			accum += f(i);
	
	return accum;
}

///////////////////////////////////////////////////////////////////////////////
// MPI functions
///////////////////////////////////////////////////////////////////////////////
void calculate_indexes (int n, int commsize)
{
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

void scatter_indexes ()
{
	// Scatter 'from' and 'to' indexes to all ranks
	MPI_Scatter (gi, 2, MPI_INT, li, 2, MPI_INT, 0, MPI_COMM_WORLD);
}

void worker ()
{
	// Calculate my range of elements
	mycalculation = integral (li[0], li[1], &zeta);
}

double reducesum (int myrank, int commsize)
{
#ifdef MPISUM
	// Unused
	(void)(myrank);
	(void)(commsize);
	
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

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////
int vtest (void)
{
	FILE *fp = fopen ("test-results.txt", "w");

	// test stuff
	
	fclose (fp);
	return 0;
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
#ifdef VTEST
	return vtest ();
#else
	// MPI Init
	int commsize;
	int myrank;
	MPI_Init (0, 0);
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
	
	// Globals allocation and index generation
	if (!myrank)
	{
		// Allocate global indexes
		gi = malloc (commsize * 2 * sizeof (int));

		// Allocate vector of elements to sum
		numbers = malloc (commsize * sizeof (double));

		// Generate indexes
		calculate_indexes (n, commsize);
	}

	// Distribute indexes
	scatter_indexes ();

	// Perform main work
	worker ();
	
	// Reduce each process process partial sum to one
	double sum = reducesum (myrank, commsize);

	// Final work for rank 0
	if (!myrank)
	{
		// Print
		printf ("%.17f\n", sqrt (sum * 6.0));

		// Free globals
		free (numbers);
		free (gi);
	}
	
	// Done
	MPI_Finalize ();
	return 0;
#endif
}
