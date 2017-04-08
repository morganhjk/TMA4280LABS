///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

///////////////////////////////////////////////////////////////////////////////
// Variables
///////////////////////////////////////////////////////////////////////////////
// These define this process working area
int from;
int to;

// These are for everyone
int *fromarray;
int *toarray;
int *adderarray;

// Stuff kept from buildlists for debugging
int list_div;
int list_rem;
int list_m;
int list_nprocs;

///////////////////////////////////////////////////////////////////////////////
// Functions
///////////////////////////////////////////////////////////////////////////////
void buildlists (int m, int nprocs, int myrank)
{
	// Allocate space
	fromarray = (int*) malloc (nprocs * sizeof (int));
	toarray = (int*) malloc (nprocs * sizeof (int));
	adderarray = (int*) malloc (nprocs * sizeof (int));

	// Calculate the number of rows per process and the remainder
	list_div = m / nprocs;
	list_rem = m % nprocs;

	// Calculate to and from row index for each process
	for (int i = 0; i < nprocs; i++)
	{
		// Calculate offset due to remainder in previous processes
		int offset = i < list_rem ? i : list_rem;

		// Calculate an adder for the number of rows we need to do
		adderarray[i] = i < list_rem ? 1 : 0;

		// Calculate the start
		fromarray[i] = i * list_div + offset;

		// Calculate the end
		toarray[i] = fromarray[i] + list_div + adderarray[i];
	}

	// Set my from and to
	from = fromarray[myrank];
	to = toarray[myrank];

	// Store m and nprocs
	list_m = m;
	list_nprocs = nprocs;
}

void destroylists ()
{
	free (fromarray);
	free (toarray);
	free (adderarray);
}

int getfromrow (int rank)
{
	return fromarray[rank];
}

int gettorow (int rank)
{
	return toarray[rank];
}

int getadder (int rank)
{
	return adderarray[rank];
}

void printlistdebug (int myrank)
{
	if (!myrank)
	{
		printf ("div %i rem %i m %i nprocs %i\n", list_div, list_rem, list_m, list_nprocs);

		for (int i = 0; i < list_nprocs; i++)
			printf ("Rank %i calculating from %i to %i\n", i, fromarray[i], toarray[i]);
	}
}

void testlist (void)
{
	// Rebuild lists
	destroylists ();
	buildlists (31, 3, 0);

	// Check all values
	if (getfromrow (0) != 0) printf ("testlist failed for getfromrow (0)\n");
	if (getfromrow (1) != 11) printf ("testlist failed for getfromrow (1)\n");
	if (getfromrow (2) != 21) printf ("testlist failed for getfromrow (2)\n");

	if (gettorow (0) != 11) printf ("testlist failed for gettorow (0)\n");
	if (gettorow (1) != 21) printf ("testlist failed for gettorow (1)\n");
	if (gettorow (2) != 31) printf ("testlist failed for gettorow (2)\n");

	if (list_div != 10) printf ("testlist failed for div\n");
	if (list_rem != 1) printf ("testlist failed for rem\n");
	if (list_m != 31) printf ("testlist failed for m\n");
	if (list_nprocs != 3) printf ("testlist failed for nprocs\n");

	// Print ok message
	printf ("Testing list done\n");

	// Free lists
	destroylists ();

	// Exit
	MPI_Finalize ();
	exit(0);
}
