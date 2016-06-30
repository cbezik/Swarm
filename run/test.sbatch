#!/bin/bash

#SBATCH --job-name=test_sbatch
#SBATCH --nodes=24
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --partition=depablo
#SBATCH --share
#SBATCH --output=mpi_test_1.dat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bezik@uchicago.edu

module load openmpi

cd /home/bezik

make clean
make

mpiexec -np 24 ./mpi_hello_world
wait
