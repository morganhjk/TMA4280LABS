///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

///////////////////////////////////////////////////////////////////////////////
// External variables
///////////////////////////////////////////////////////////////////////////////
extern int commsize;
extern int myrank;

///////////////////////////////////////////////////////////////////////////////
// Variables
///////////////////////////////////////////////////////////////////////////////
// These define this process working area
int from;
int to;

// These are for everyone
int *fromarray;
int *toarray;

///////////////////////////////////////////////////////////////////////////////
// Functions
///////////////////////////////////////////////////////////////////////////////
void definearea (int m)
{
	// Calculate the number of rows per process and the remainder
	int div = m / commsize;
	int rem = m % commsize;

	// Calculate offset due to remainder in previous processes
	int offset = myrank < rem ? myrank * 1 : rem * 1;

	// Calculate an adder for the number of rows we need to do
	int adder = myrank < rem ? 1 : 0;

	// Calculate the start
	from = myrank * div + offset;

	// Calculate the end
	to = from + div + adder;

	// Debug
#if 0
	if (!myrank)
		printf ("div %i rem %i m %i commsize %i\n", div, rem, m, commsize);

	printf ("Process %i calculating from %i to %i offset %i adder %i\n", myrank, from, to, offset, adder);
#endif
}

void buildlists (int m, int nprocs)
{
	// Allocate space
	fromarray = (int*) malloc (nprocs * sizeof (int));
	toarray = (int*) malloc (nprocs * sizeof (int));
	
	// Calculate the number of rows per process and the remainder
	int div = m / nprocs;
	int rem = m % nprocs;
	
	// Calculate to and from row index for each process
	for (int i = 0; i < nprocs; i++)
	{
		// Calculate offset due to remainder in previous processes
		int offset = i < rem ? i : rem;

		// Calculate an adder for the number of rows we need to do
		int adder = i < rem ? 1 : 0;

		// Calculate the start
		fromarray[i] = i * div + offset;

		// Calculate the end
		toarray[i] = fromarray[i] + div + adder;
	}
	
	// Debug
#if 0
	if (!myrank)
	{
		printf ("div %i rem %i m %i nprocs %i\n", div, rem, m, nprocs);
		
		for (int i = 0; i < nprocs; i++)
			printf ("Rank %i calculating from %i to %i\n", i, fromarray[i], toarray[i]);
	}
#endif
}

void destroylists ()
{
	free (fromarray);
	free (toarray);
}

int getfromrow (int rank)
{
	return fromarray[rank];
}

int gettorow (int rank)
{
	return toarray[rank];
}
