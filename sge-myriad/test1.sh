#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=72:0:0
#$ -l mem=2G
#$ -l tmpfs=4G
#$ -pe mpi 36
#$ -N test1
#$ -cwd
#$ -o output/2021.Test1/test1.out

cd /home/uctpln0/FruitDemand/code/fortran/source
hostname
date
mpirun -genv I_MPI_FABRICS=tcp -np $NSLOTS ./SparseDemand_mpi.exe ../inputs/Test1/A1a.estimate.prop
date
