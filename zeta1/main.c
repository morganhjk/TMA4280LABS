#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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

int main (int argc, char **argv)
{
	// MPI Init
	int commsize;
	int myrank;
	MPI_Init (0, 0);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	// Check input arguments
	if (argc != 2 || ((commsize & (commsize - 1)) != 0))
	{
		if (!myrank)
			printf ("usage: ./zeta1 n, where n is a power of two\n");
		return 0;
	}
	
	int n = atoi (argv[1]);
	
	if (n <= 0)
	{
		if (!myrank)
			printf ("n is bullshit, try again with different n\n");
		return 2;
	}
	
	// Work some
	double ret = 0.0;
	int ndiv = n / commsize;
	int nrem = n % commsize;
	int nsize = commsize + (nrem ? 1 : 0);
	
	if (!myrank)
	{
		double *numbers = malloc (nsize * sizeof (double));
		
		numbers[0] = integral (1, ndiv, &zeta);
		
		if (nrem)
			numbers[nsize-1] = integral (ndiv * commsize + 1, ndiv * commsize + nrem, &zeta);
		
		for (int i = 1; i < commsize; i++)
		{
			MPI_Status s;
			MPI_Recv (&numbers[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &s);
		}
		
		for (int i = 0; i < nsize; i++)
			ret += numbers[i];
		
		free (numbers);
	}
	else
	{
		ret = integral (ndiv * myrank + 1, ndiv * (myrank+1), &zeta);
		MPI_Send (&ret, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
	}
	
	// Print
	if (!myrank)
		printf ("%.17f\n", sqrt (ret * 6.0));
	
	// Done
	MPI_Finalize ();
	return 0;
}
