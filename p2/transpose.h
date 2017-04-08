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
 * Transpose matrix b in parallel, and store it in bt.
 */
void transpose (double **bt, double **b, size_t m);

/*
 * For debug purposes:
 * Print the contents of the arrays made by inittranspose.
 */
void printtranspose (int rank, int m);

/*
 * For debug purposes:
 * Print matrix of size m in c.
 */
void printtestmatrix (int rank, double *c, int m);

/*
 * Test the transpose stuff. Note that it will destroy and rebuild
 * lists as it sees fit and may not be used in production.
 *
 * This test hast to be run with 3 processes. It allocates an m=7
 * matrix, initializes it with values unique to each process and
 * tries to distribute/transpose it across all three processes.
 * Each process then prints what it got to standard output for
 * visual inspection.
 */
void testtranspose (void);

#endif
