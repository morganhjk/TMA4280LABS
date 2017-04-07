/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 *
 * Parallel version by Morgan Kristiansen and Erlend Sveen, 2017
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
// Definitions
///////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

#define PRINTTABLES 0

///////////////////////////////////////////////////////////////////////////////
// Function prototypes
///////////////////////////////////////////////////////////////////////////////
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);

///////////////////////////////////////////////////////////////////////////////
// Fortran functions
///////////////////////////////////////////////////////////////////////////////
// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

///////////////////////////////////////////////////////////////////////////////
// Global variables
///////////////////////////////////////////////////////////////////////////////
// Size and rank
int commsize;
int myrank;

// Number of threads per process
int nthreads = 1;

/*
 *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
 *  conditions are applied on the boundary:
 *  - the number of grid points in each direction is n+1,
 *  - the number of degrees of freedom in each direction is m = n-1,
 *  - the mesh size is constant h = 1/n.
 */
int n, m;
real h;

///////////////////////////////////////////////////////////////////////////////
// Functions
///////////////////////////////////////////////////////////////////////////////
void usage (char *appname)
{
	printf ("usage: mpirun -np <p> %s <n> <t>\n", appname);
	printf ("where:\n");
	printf (" <p> is a positive integer (number of processes)\n");
	printf (" <n> is a power of two (problem size)\n");
	printf (" <t> is a positive integer (number of threads per process)\n");
}

void initialize (int argc, char **argv)
{
	// MPI Init
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Check number of arguments
	if (argc != 3)
	{
		if (!myrank)
			usage (argv[0]);
		exit (0);
	}

	// Get problem size (n)
	n = atoi (argv[1]);

	// Get number of threads (nthreads)
	nthreads = atoi (argv[2]);

	// Check problem size
	if (((n & (n - 1)) != 0))
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a power of two, try again with a different <p>\n", n);
			usage (argv[0]);
		}
		exit (1);
	}

	// Check n
	if (n <= 0)
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a positive integer, try again with a different <n>\n", n);
			usage (argv[0]);
		}
		exit (2);
	}

	// Check nthreads
	if (nthreads <= 0)
	{
		if (!myrank)
		{
			fprintf (stderr, "'%i' is not a positive integer, try again with a different <t>\n", nthreads);
			usage (argv[0]);
		}
		exit (3);
	}

	// Set number of degrees of freedom and mesh size
	m = n - 1;
	h = 1.0 / n;
}

void findmax (real **b, int m)
{
	// Global max and my max
	double u_max = 0.0;
	double mymax = 0.0;

	// Calculate
	for (size_t i = from; i < to; i++)
		for (size_t j = 0; j < m; j++)
			mymax = mymax > b[i][j] ? mymax : b[i][j];

	// Reduce max
	MPI_Allreduce (&mymax, &u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	// Print max if we're process 0
	if (!myrank)
		printf ("u_max = %e\n", u_max);
}

void printresults (real **bt, real **b, int m)
{
#if PRINTTABLES
	// Print the matrix from rank zero
	if (!myrank)
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
				printf (" %f", b[i][j]);

			printf ("\n");
		}
	}
#endif
}

real checkfunc (real x, real y)
{
	return sin (PI * x) * sin (2.0 * PI * y);
}

void errorcheck (real **b, real *grid, int m)
{
#if PRINTTABLES
	if (!myrank)
		printf ("\n");
#endif

	real **u = mk_2D_array (m, m, true);
	
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < m; j++)
			u[i][j] = checkfunc (grid[i+1], grid[j+1]);
	
	double maxerr = 0.0;
	
	if (!myrank)
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				double error = u[i][j] - b[i][j];
				
#if PRINTTABLES
				printf (" %f", error);
#endif
				
				if (error < 0.0)
					error = -error;
				
				if (error > maxerr)
					maxerr = error;
			}

#if PRINTTABLES
			printf ("\n");
#endif
		}
		
		printf ("Largest error encountered: %f\n", maxerr);
	}
	
	free (u);
}

int main (int argc, char **argv)
{
	// Prepare ourselves for parallel operation
	initialize (argc, argv);

	// Define this process working area
	buildlists (m, commsize, myrank);
	inittranspose (myrank, m);

	/*
	 * Grid points are generated with constant mesh size on both x- and y-axis.
	 */
	real *grid = mk_1D_array(n+1, false);

#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = 0; i < n+1; i++)
		grid[i] = i * h;

	/*
	 * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
	 * defined Chapter 9. page 93 of the Lecture Notes.
	 * Note that the indexing starts from zero here, thus i+1.
	 */
	real *diag = mk_1D_array(m, false);

#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = 0; i < m; i++)
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));

	/*
	 * Allocate the matrices b and bt which will be used for storing value of
	 * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
	 */
	real **b = mk_2D_array(m, m, false);
	real **bt = mk_2D_array(m, m, false);

	/*
	 * This vector will holds coefficients of the Discrete Sine Transform (DST)
	 * but also of the Fast Fourier Transform used in the FORTRAN code.
	 * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
	 * - Fourier coefficients are complex so storage is used for the real part
	 *   and the imaginary part.
	 * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
	 *   DST coefficients are defined for j [[ 0, n-1 ]].
	 * As explained in the Lecture notes coefficients for positive j are stored
	 * first.
	 * The array is allocated once and passed as arguments to avoid doings 
	 * reallocations at each function call.
	 */
	int nn = 4 * n;
	real *z[nthreads];

	// Create a vector for each thread in this process
	for (int i = 0; i < nthreads; i++)
		z[i] = mk_1D_array(nn, false);

	/*
	 * Initialize the right hand side data for a given rhs function.
	 * Note that the right hand-side is set at nodes corresponding to degrees
	 * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
	 * 
	 */
#pragma omp parallel for num_threads(nthreads) schedule(static) collapse(2)
	for (size_t i = from; i < to; i++)
		for (size_t j = 0; j < m; j++)
			b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);

	/*
	 * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
	 * Instead of using two matrix-matrix products the Discrete Sine Transform
	 * (DST) is used.
	 * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
	 * The array zz is used as storage for DST coefficients and internally for 
	 * FFT coefficients in fst_ and fstinv_.
	 * In functions fst_ and fst_inv_ coefficients are written back to the input 
	 * array (first argument) so that the initial values are overwritten.
	 */
#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = from; i < to; i++)
		fst_ (b[i], &n, z[omp_get_thread_num()], &nn);

	transpose (bt, b, m);

#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = from; i < to; i++)
		fstinv_ (bt[i], &n, z[omp_get_thread_num()], &nn);

	// Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
#pragma omp parallel for num_threads(nthreads) schedule(static) collapse(2)
	for (size_t i = from; i < to; i++)
		for (size_t j = 0; j < m; j++)
			bt[i][j] = bt[i][j] / (diag[i] + diag[j]);

	// Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = from; i < to; i++)
		fst_ (bt[i], &n, z[omp_get_thread_num()], &nn);

	transpose (b, bt, m);

#pragma omp parallel for num_threads(nthreads) schedule(static)
	for (size_t i = from; i < to; i++)
		fstinv_ (b[i], &n, z[omp_get_thread_num()], &nn);

	// Compute maximal value of solution for convergence analysis in L_\infty norm.
	findmax (b, m);

	// Cheap hack: transpose twice to sync. Following functions are not parallel.
	transpose (bt, b, m);
	transpose (b, bt, m);

	// Print result matrix
	printresults (bt, b, m);

	// Error check as in appendix B
	errorcheck (b, grid, m);

	// Done
	deinittranspose ();
	destroylists ();
	MPI_Finalize ();
	return 0;
}

/*
* This function is used for initializing the right-hand side of the equation.
* Other functions can be defined to swtich between problem definitions.
*/

real rhs (real x, real y)
{
#if 0
	return 2 * (y - y*y + x - x*x);
#else
	return 5.0 * PI * PI * sin (PI * x) * sin (2.0 * PI * y);
#endif
}

/*
* The allocation of a vectore of size n is done with just allocating an array.
* The only thing to notice here is the use of calloc to zero the array.
*/

real *mk_1D_array (size_t n, bool zero)
{
	if (zero)
		return (real*) calloc (n, sizeof(real));

	return (real*) malloc (n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array (size_t n1, size_t n2, bool zero)
{
	// 1
	real **ret = (real**) malloc (n1 * sizeof(real*));

	// 2
	if (zero)
		ret[0] = (real*) calloc (n1 * n2, sizeof(real));
	else
		ret[0] = (real*) malloc (n1 * n2 * sizeof(real));

	// 3
	for (size_t i = 1; i < n1; i++)
		ret[i] = ret[i-1] + n2;

	return ret;
}
