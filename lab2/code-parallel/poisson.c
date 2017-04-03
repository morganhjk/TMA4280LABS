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

///////////////////////////////////////////////////////////////////////////////
// Definitions
///////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

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

// These define this process working area
int from;
int to;

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
	int n = atoi (argv[1]);

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
}

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
#if 1
	if (!myrank)
		printf ("div %i rem %i m %i commsize %i\n", div, rem, m, commsize);

	printf ("Process %i calculating from %i to %i offset %i adder %i\n", myrank, from, to, offset, adder);
#endif
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

int main (int argc, char **argv)
{
	// Prepare ourselves for parallel operation
	initialize (argc, argv);

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    int n = atoi(argv[1]);
    int m = n - 1;
    real h = 1.0 / n;

	// Define this process working area
	definearea (m);

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);

    for (size_t i = 0; i < n+1; i++)
        grid[i] = i * h;

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);

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
    real *z = mk_1D_array(nn, false);

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     * 
     */
    for (size_t i = 0; i < m; i++)
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
    for (size_t i = 0; i < m; i++)
        fst_(b[i], &n, z, &nn);

    transpose(bt, b, m);

    for (size_t i = 0; i < m; i++)
        fstinv_(bt[i], &n, z, &nn);

	// Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);

	// Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
    for (size_t i = 0; i < m; i++)
        fst_(bt[i], &n, z, &nn);

    transpose(b, bt, m);

    for (size_t i = 0; i < m; i++)
        fstinv_(b[i], &n, z, &nn);

    // Compute maximal value of solution for convergence analysis in L_\infty norm.
    findmax (b, m);

	// Done
	MPI_Finalize ();
	return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs (real x, real y)
{
    return 2 * (y - y*y + x - x*x);
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void transpose (real **bt, real **b, size_t m)
{
#if 0
	// Allocate stuff
	int *sendcounts = (int*) malloc (commsize * sizeof (int));
	int *recvcounts = (int*) malloc (commsize * sizeof (int));
	int *senddispl  = (int*) malloc (commsize * sizeof (int));
	int *recvdispl  = (int*) malloc (commsize * sizeof (int));
	
	// Init
	for (int i = 0; i < commsize; i++)
	{
		// At this point, from and to describe the rows for this process
		// Get n
		int nrows = to - from;
		
		sendcounts[i] = 1;
		recvcounts[i] = 1;
		senddispl[i] = 1;
		recvdispl[i] = 1;
	}
	
	// Send
	MPI_Alltoall (b[0], m/commsize, MPI_DOUBLE,
		bt[0], m/commsize, MPI_DOUBLE, MPI_COMM_WORLD);
	
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
#else
	// Naive
	for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < m; j++)
            bt[i][j] = b[j][i];
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
