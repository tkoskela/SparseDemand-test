#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/IO/Lewbel/code/fortran/source
#$ -o /home/uctpln0/IO/Lewbel/code/fortran/output/a1_mpi.out
#$ -N a1_mpi
#$ -j y
#$ -cwd
#$ -l tmem=2G,h_vmem=2G
#$ -l h_rt=8:0:0
#$ -pe orte 16

cd /home/uctpln0/IO/Lewbel/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
date
#./SparseDemand.dbg ../inputs/A1.prop
mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A1.prop
date
