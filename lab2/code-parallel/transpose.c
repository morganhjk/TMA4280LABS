

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "list.h"

extern int commsize;
extern int myrank;
extern int from;
extern int to;

int *sendcounts;
int *recvcounts;
int *senddispl;
int *recvdispl;

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
		senddispl[i] = from * m;
		
		// Get number of rows to receive
		int nrowsrecv = gettorow(i) - getfromrow(i);
		
		// Set amount and receive index
		recvcounts[i] = nrowsrecv * m;
		recvdispl[i] = getfromrow(i) * m;
	}
}

void printtranspose (int rank, int m)
{
	// Print debug info
	for (int i = 0; i < commsize; i++)
	{
		int nrows = gettorow(rank) - getfromrow(rank);
		printf ("from %i to rank %i nrows %i sendcounts[i] %i senddispl[i] %i recvcounts[i] %i recvdispl[i] %i\n",
			rank, i, nrows, sendcounts[i], senddispl[i], recvcounts[i], recvdispl[i]);
		
	}
}

void transpose (double **bt, double **b, size_t m)
{
#if 1
	// Init
	inittranspose (myrank, m);
	printtranspose (myrank, m);
	
	// Send
	//MPI_Alltoall (b[0], m/commsize, MPI_DOUBLE,
	//	bt[0], m/commsize, MPI_DOUBLE, MPI_COMM_WORLD);
	
	/*void parallel_transpose_all(int id,double *lphi,double *tlphi,double *tempaux1, int m,int ln){
        serial_transpose(lphi, tempaux1,ln,m);
        MPI_Datatype    block_t;
        MPI_Datatype    ub_block_t;   //upperbound row
        int             bl[2];      //block length
        MPI_Aint        dl[2];      //displacement location
        MPI_Datatype    type[2]; 
        MPI_Type_vector(ln, ln, m, MPI_DOUBLE, &block_t);
        bl[0] = bl[1] = 1;
        dl[0] = 0;
        dl[1] = ln*sizeof(double);
        type[0] = block_t;
        type[1] = MPI_UB;
        MPI_Type_struct(2, bl, dl, type, &ub_block_t);
        MPI_Type_commit(&ub_block_t);
 
       MPI_Alltoall( tempaux1,1,ub_block_t,tlphi,(ln*ln),MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Type_free(&ub_block_t);
        MPI_Type_free(&block_t);
}*/

	// Free
	free (sendcounts);
	free (recvcounts);
	free (senddispl);
	free (recvdispl);
	
	MPI_Finalize ();
	exit(0);
#else
	// Naive
	for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            bt[i][j] = b[j][i];
#endif
}

