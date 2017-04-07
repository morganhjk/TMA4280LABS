/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

///////////////////////////////////////////////////////////////////////////////
// Includes
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "list.h"
#include "transpose.h"

///////////////////////////////////////////////////////////////////////////////
// Extern variables
///////////////////////////////////////////////////////////////////////////////
extern int commsize;
extern int myrank;

///////////////////////////////////////////////////////////////////////////////
// Private variables
///////////////////////////////////////////////////////////////////////////////
// Send counts and displacements
int *sendcounts;
int *recvcounts;
int *senddispl;
int *recvdispl;

// MPI Datatypes
MPI_Datatype matrixcolumn;
MPI_Datatype matrixcolumntype;

// MPI Datatypes for alternative transpose method
MPI_Datatype sendtypelarge;
MPI_Datatype sendtypesmall;
MPI_Datatype recvtypelarge;
MPI_Datatype recvtypesmall;

///////////////////////////////////////////////////////////////////////////////
// Functions
///////////////////////////////////////////////////////////////////////////////
void inittranspose (int rank, int m)
{
	// Allocate stuff
	sendcounts = (int*) malloc (commsize * sizeof (int));
	recvcounts = (int*) malloc (commsize * sizeof (int));
	senddispl  = (int*) malloc (commsize * sizeof (int));
	recvdispl  = (int*) malloc (commsize * sizeof (int));
	
	// Init each entry
	for (int i = 0; i < commsize; i++)
	{
		// Get the number of rows to send
		int nrows = gettorow(rank) - getfromrow(rank);
		
		// Set our amount and displacement. We send all our data to everyone
		sendcounts[i] = nrows * m;
		senddispl[i] = getfromrow(rank) * m;
		
		// Get number of rows to receive
		int nrowsrecv = gettorow(i) - getfromrow(i);
		
		// Set amount and receive index
		recvcounts[i] = nrowsrecv;
		recvdispl[i] = getfromrow(i);
	}
	
	// Create MPI data type for matrix columns
	MPI_Type_vector (m, 1, m, MPI_DOUBLE, &matrixcolumn);
	MPI_Type_commit (&matrixcolumn);
	MPI_Type_create_resized (matrixcolumn, 0, sizeof(double), &matrixcolumntype);
	MPI_Type_commit (&matrixcolumntype);
	
#if 0
	//
	int rowstosend = gettorow(rank) - getfromrow(rank);
	int ncolslarge = gettorow(0) - getfromrow(0);
	int ncolssmall = gettorow(commsize-1) - getfromrow(commsize-1);

	//
	MPI_Type_vector (rowstosend, ncolslarge, m, MPI_DOUBLE, &sendtypelarge);
	MPI_Type_vector (rowstosend, ncolssmall, m, MPI_DOUBLE, &sendtypesmall);
	MPI_Type_commit (&sendtypelarge);
	MPI_Type_commit (&sendtypesmall);

	//
	MPI_Type_vector (rowstosend, ncolslarge, m, MPI_DOUBLE, &recvtypelarge);
	MPI_Type_vector (rowstosend, ncolssmall, m, MPI_DOUBLE, &recvtypesmall);
	MPI_Type_commit (&recvtypelarge);
	MPI_Type_commit (&recvtypesmall);
#endif
}

void deinittranspose (void)
{
	free (sendcounts);
	free (recvcounts);
	free (senddispl);
	free (recvdispl);
}

void transpose (double **bt, double **b, size_t m)
{
#if 1
	// Parallel version
	MPI_Alltoallv (b[0], sendcounts, senddispl, MPI_DOUBLE,
		bt[0], recvcounts, recvdispl, matrixcolumntype, MPI_COMM_WORLD);
#endif

#if 0
	// Parallel version, timed
	double start = MPI_Wtime ();

	MPI_Alltoallv (b[0], sendcounts, senddispl, MPI_DOUBLE,
		bt[0], recvcounts, recvdispl, matrixcolumntype, MPI_COMM_WORLD);

	double end = MPI_Wtime ();

	if (!myrank)
		printf ("Transpose time %f seconds\n", (end - start));
#endif

#if 0
	// This is just a test, do not use
	double start = MPI_Wtime ();

	for (int i = 0; i < commsize; i++)
	{
		int indexto = (myrank - i) % commsize;
		
		int offset = getfromrow (myrank) * m;
		
		offset += indexto;
		
		MPI_Status status;
		
		MPI_Sendrecv (&b[0][offset], 1, (getadder(i) ? sendtypelarge : sendtypesmall),
			i, 0,
			&b[0][offset], 1, (getadder(i) ? recvtypelarge : recvtypesmall),
			i, 0, MPI_COMM_WORLD, &status);
	}

	double end = MPI_Wtime ();

	if (!myrank)
		printf ("Transpose time %f seconds\n", (end - start));
#endif

#if 0
	// Naive approach from the naive version
	for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            bt[i][j] = b[j][i];
#endif
}


///////////////////////////////////////////////////////////////////////////////
// Debug Functions
///////////////////////////////////////////////////////////////////////////////
void printtranspose (int rank, int m)
{
	// Print debug info
	for (int i = 0; i < commsize; i++)
	{
		// Recalculate nrows as it is not stored explicitly
		int nrows = gettorow(rank) - getfromrow(rank);
		
		// Print contents
		printf ("from %i to rank %i nrows %i"
				" sendcounts[i] %i senddispl[i] %i recvcounts[i] %i recvdispl[i] %i\n",
			rank, i, nrows, sendcounts[i], senddispl[i], recvcounts[i], recvdispl[i]);
		
	}
}

void printtestmatrix (int rank, double *c, int m)
{
	if (rank != myrank)
		return;

	if (m == 3)
	printf ("rank %i "
			"\t%f %f %f\n"
			"\t%f %f %f\n"
			"\t%f %f %f\n",
			rank,
			c[0], c[1], c[2],
			c[3], c[4], c[5],
			c[6], c[7], c[8]);
	else
	printf ("rank %i "
			"\t%f %f %f %f %f %f %f\n"
			"\t%f %f %f %f %f %f %f\n"
			"\t%f %f %f %f %f %f %f\n",
			rank,
			c[0], c[1], c[2], c[3], c[4], c[5], c[6],
			c[7], c[8], c[9], c[10], c[11], c[12], c[13],
			c[14], c[15], c[16], c[17], c[18], c[19], c[20]);
}

void testtranspose ()
{
	// Check that we are running with correct number of processes first
	if (commsize != 3)
	{
		if (!myrank)
			printf ("incorrect commsize, must test with 3\n");
		MPI_Finalize ();
		exit (0);
	}
	
	// Test matrix size
	int m = 7;
	
	// Rebuild lists for our test purposes
	destroylists ();
	buildlists (m, commsize);
	
	// Create a source and destination matrix
	double *b = (double*) malloc (m * m * sizeof (double));
	double *c = (double*) malloc (m * m * sizeof (double));
	
	// Initialize source matrix with something we recognize, and set destination to zero
	for (int i = 0; i < m*m; i++)
	{
		b[i] = (myrank * 100) + i;
		c[i] = 0;
	}
	
	// Initialize Alltoallv arrays
	inittranspose (myrank, m);
	
	// Print their contents for examination
	printtranspose (myrank, m);
	
	// Do the actual transmission and transposition using a single mpi call
	transpose (&c, &b, m);
	
	// Print the results from each process
	for (int i = 0; i < commsize; i++)
	{
		MPI_Barrier (MPI_COMM_WORLD);
		printtestmatrix (i, c, m);
	}
	
	// Free transpose arrays
	deinittranspose ();
	
	// Free matrix arrays
	free (c);
	free (b);
	
	// Free lists
	destroylists ();
	
	// Exit
	MPI_Finalize ();
	exit(0);
}
