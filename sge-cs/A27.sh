#$ -S /bin/bash
#$ -wd /home/uctpln0/FuirtDemand/code/fortran/sge_cs
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A27.out
#$ -N fruit27
#$ -j y
#$ -cwd
#$ -l tmem=4G,h_vmem=4G,hostname=burns*
#$ -l h_rt=240:0:0
## ,hostname=burns-*
##$ -l hostname=burns*\|zeppo*\|fry*\|larry*\|cheech*
##$ -l hostname=burns*\|fry*\
#$ -R y
#$ -pe orte 100

cd /home/uctpln0/FruitDemand/code/fortran/source
hostname
date
which mpirun
source
mpirun -perhost 1 -np $NSLOTS -env I_MPI_FABRICS tcp ./SparseDemand_mpi.exe ../inputs/A27_new/A27A.prop
date

