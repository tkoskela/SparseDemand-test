#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/IO/Lewbel/code/fortran/source
#$ -o /home/uctpln0/IO/Lewbel/code/fortran/output/a31_5.out
#$ -N a31_5
#$ -j y
#$ -cwd
#$ -l tmem=2G,h_vmem=2G
#$ -l h_rt=12:0:0
#$ -pe orte 4

cd /home/uctpln0/IO/Lewbel/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
date
#./SparseDemand.dbg ../inputs/A1.prop
#mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A31_5.prop
mpirun -np $NSLOTS SparseDemand_mpi.dbg ../inputs/A31_5.prop
date
