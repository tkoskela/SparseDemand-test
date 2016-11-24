#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/FuirtDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A24A_mpi.out
#$ -N A24A_MPI
#$ -j y
#$ -cwd
#$ -l tmem=2G,h_vmem=2G
#$ -l h_rt=24:0:0,hostname=burns-*
#$ -pe orte 40

cd /home/uctpln0/FruitDemand/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
. /share/apps/econ/Modules/modules.sh
module load gcc/5.2.0 nag/fll6i25dc intel/composer/2015.1.133 openmpi/intel/1.10.0 totalview knitro
date
#./SparseDemand.dbg ../inputs/A1.prop
#mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A31_5.prop
mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A24/A24A.prop
date
