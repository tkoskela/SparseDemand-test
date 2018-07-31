#$ -S /bin/bash
#$ -wd /home/uctpln0/FruitDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A27A_analyse.out
#$ -N A27A_anal
#$ -j y
#$ -cwd
#$ -l tmem=4G,h_vmem=4G
#$ -l h_rt=240:0:0
##$ -l hostname=burns*,zeppo*,fry*,larry*,cheech*
#$ -R y
##$ -pe orte 201

cd /home/uctpln0/FruitDemand/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
. /share/apps/econ/Modules/modules.sh
module load gcc/5.2.0 nag/fll6i25dc intel/composer/2015.1.133 openmpi/intel/1.10.0 totalview knitro
hostname
date
which mpirun
mpirun --mca btl_tcp_if_include 10.0.0.0/8 -np $NSLOTS SparseDemand_mpi_new.exe ../inputs/A27/A27A.prop
date

