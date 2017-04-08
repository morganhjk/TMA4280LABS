#!/bin/bash
#PBS -N poisson_p26-t01
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=20:mpiprocs=13
#PBS -m abe
#PBS -M morganhk@stud.ntnu.no
cd "$PBS_O_WORKDIR" || exit 1
module load gcc openmpi

# Define executable
poisson="../build/poisson"

# Define problem sizes
nsizes=( 2048 4096 8192 16384 )

# Run program for all sizes
for nsize in "${nsizes[@]}"
do
	echo "=========================="
	time mpirun "$poisson" "$nsize" 1
	echo "=========================="
done
