#!/bin/bash
#
# A script to submit SGE_SparseDemand.exe
# Run using 'qsub SGE_SparseDemand.sh

#
# Specify the Working Directory
#$ -wd /home/uctpln0/IO/Lewbel/code/fortran
 
# Specify an output file
#$ -o /home/uctpln0/IO/Lewbel/code/fortran/source/output/a1.out

# Specify whether or not the standard error stream of the job i
# is merged into the standard output stream
# y = yes
# n = no
# $ -j y

# Specify the shell to use
#$ -S /bin/bash

# Specify a job name for the job
#$ -N SparseDemand

# Launch the job in a queue meeting the given resource request list
#$ -l tmem=1G,h_vmem=1G
#$ -l h_rt=8:0:0

# Specify a parallel environment
#$ -pe orte 16
#$ -V

source /home/uctpln0/GeneralCode/config/SetEnv.sh

cd /home/uctpln0/IO/Lewbel/code/fortran/source

date
mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A1.prop
date
