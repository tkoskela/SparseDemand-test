#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=72:0:0
#$ -l mem=2G
#$ -l tmpfs=4G
#$ -pe mpi 400
#$ -N fruit1
#$ -cwd
#$ -o output/2021.04/fruit1.out

cd /home/uctpln0/FruitDemand/code/fortran/source
hostname
date
mpirun -genv I_MPI_FABRICS=tcp -np $NSLOTS ./SparseDemand_mpi.exe ../inputs/A27.2021.04/A27A.prop
#mpirun -genv I_MPI_FABRICS=tcp -np $NSLOTS ./SparseDemand_mpi.dbg ../inputs/A27.2021.04/A27A.prop
date
