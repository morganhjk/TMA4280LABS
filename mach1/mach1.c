#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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
	double reta = 0.0;
	double retb = 0.0;
	int ndiv = n / commsize;
	int nrem = n % commsize;
	int nsize = commsize + (nrem ? 1 : 0);
	
	if (!myrank)
	{
		double *numbersa = malloc (nsize * sizeof (double));
		double *numbersb = malloc (nsize * sizeof (double));
		
		numbersa[0] = integral (1, ndiv, (1.0 / 5.0), &machin);
		numbersb[0] = integral (1, ndiv, (1.0 / 239.0), &machin);
		
		if (nrem)
		{
			numbersa[nsize-1] = integral (ndiv * commsize + 1, ndiv * commsize + nrem, (1.0 / 5.0), &machin);
			numbersb[nsize-1] = integral (ndiv * commsize + 1, ndiv * commsize + nrem, (1.0 / 239.0), &machin);
		}
		
		for (int i = 1; i < commsize; i++)
		{
			MPI_Status s;
			MPI_Recv (&numbersa[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &s);
			MPI_Recv (&numbersb[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &s);
		}
		
		for (int i = 0; i < nsize; i++)
		{
			reta += numbersa[i];
			retb += numbersb[i];
		}
		
		free (numbersa);
		free (numbersb);
	}
	else
	{
		reta = integral (ndiv * myrank + 1, ndiv * (myrank+1), (1.0 / 5.0), &machin);
		retb = integral (ndiv * myrank + 1, ndiv * (myrank+1), (1.0 / 239.0), &machin);
		MPI_Send (&reta, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
		MPI_Send (&retb, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
	}
	
	// Print
	if (!myrank)
		printf ("%.17f\n", 4.0 * (4 * reta - retb));
	
	// Done
	MPI_Finalize ();
	return 0;
}
