#ifndef TRANSPOSE_H
#define TRANSPOSE_H

/*
 * Before we can transpose anything, we must initialize the arrays
 * and MPI datatypes that MPI_Alltoallv needs. Takes the current
 * process' rank and matrix dimension.
 */
void inittranspose (int rank, int m);

/*
 * Free arrays made by inittranspose.
 */
void deinittranspose (void);

/*
 * For debug purposes: Print the contents of the arrays made by
 * inittranspose.
 */
void printtranspose (int rank, int m);

/*
 * Print matrix of fixed size 3x3 in c, if rank equals the current
 * process' rank.
 */
void printtestmatrix (int rank, double *c);

/*
 * Test the transpose stuff. Note that it will destroy and rebuild
 * lists as it sees fit and may not be used in production.
 */
void testtranspose (void);

/*
 * Transpose matrix b in parallel, and store it in bt.
 */
void transpose (double **bt, double **b, size_t m);

#endif
