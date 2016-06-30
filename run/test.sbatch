#!/bin/bash

#SBATCH --job-name=swarm_test
#SBATCH --nodes=4
#SBATCH --time=00:10:00
#SBATCH --partition=sandyb
#SBATCH --qos=debug
#SBATCH --share
#SBATCH --output=swarm_test.dat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bezik@uchicago.edu

module unload intelmpi
module load openmpi
module load gcc/4.8

cd /home/bezik/Swarm/src

make clean
make

mpiexec -np 4 ./swarm

wait
