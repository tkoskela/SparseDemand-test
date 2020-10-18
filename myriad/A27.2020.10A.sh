#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=120:0:0
#$ -l mem=4G
#$ -l tmpfs=4G
#$ -pe mpi 200
#$ -t 1
#$ -N fruit
#$ -wd /home/uctpln0/FruitDemand/code/fortran/source
#$ -o /home/uctpln0/Scratch/FruitDemand/output/a27.2020.10.out

source /shared/ucl/apps/bin/defmods
module load nag/fortran/mark26/intel-2017

cd /home/uctpln0/FruitDemand/code/fortran/source
date
export I_MPI_FABRICS=shm,tcp
mpirun -genv I_MPI_FABRICS=tcp -np $NSLOTS ./SparseDemand_mpi.exe ../inputs/A27.2020.10/A27A.prop
date

