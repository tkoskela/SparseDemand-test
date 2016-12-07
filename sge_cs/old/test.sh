#!/bin/bash
#$ -S /bin/bash
#$ -wd /home/uctpln0/FruitDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A24B_mpi.out
#$ -N A24B_MPI
#$ -j y
#$ -R y
#$ -cwd
#$ -V
#$ -t 1
#$ -l tmem=2G,h_vmem=2G
#$ -l h_rt=240:0:0
##,hostname=burns-*
##$ -pe orte 40
##  -terse  -R y -t 1-$np -V $res_list 

cd /home/uctpln0/FruitDemand/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
. /share/apps/econ/Modules/modules.sh
module load gcc/5.2.0 nag/fll6i25dc intel/composer/2015.1.133 openmpi/intel/1.10.0 totalview knitro
date
./sleep.exe
date
