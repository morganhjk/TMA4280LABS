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
		for (int i = from; i <= to; i++)
			accum += f(i, x);
	
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
	// Declare and init local sums
	double reta = 0.0;
	double retb = 0.0;

	// Calculate my range of elements
	reta = integral (li[0], li[1], (1.0 / 5.0), &machin);
	retb = integral (li[0], li[1], (1.0 / 239.0), &machin);

	// Gather sums from all ranks
	MPI_Gather (&reta, 1, MPI_DOUBLE, numbersa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather (&retb, 1, MPI_DOUBLE, numbersb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

		// Allocate vectors of elements to sum
		numbersa = malloc (commsize * sizeof (double));
		numbersb = malloc (commsize * sizeof (double));

		// Generate indexes
		calculate_indexes (n, commsize);
	}

	// Distribute indexes
	scatter_indexes ();

	// Perform main work
	worker ();

	// Final work for rank 0
	if (!myrank)
	{
		// Sum vector elements
		double suma = 0.0;
		double sumb = 0.0;
		for (int i = 0; i < commsize; i++)
		{
			suma += numbersa[i];
			sumb += numbersb[i];
		}

		// Print
		printf ("%.17f\n", 4.0 * (4.0 * suma - sumb));

		// Free globals
		free (numbersa);
		free (numbersb);
		free (gi);
	}
	
	// Done
	MPI_Finalize ();
	return 0;
#endif
}
