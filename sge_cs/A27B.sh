#$ -S /bin/bash
#$ -wd /home/uctpln0/FuirtDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A27B.out
#$ -N fruit27B
#$ -j y
#$ -cwd
#$ -l tmem=4G,h_vmem=4G
#$ -l h_rt=240:0:0
## ,hostname=burns-*
##$ -l hostname=burns*\|zeppo*\|fry*\|larry*\|cheech*
##$ -l hostname=burns*\|fry*\
#$ -R y
#$ -pe orte 200

cd /home/uctpln0/FruitDemand/code/fortran/source
hostname
date
which mpirun
mpirun -genv I_MPI_FABRICS=tcp -np $NSLOTS ./SparseDemand_mpi.exe ../inputs/A27_new/A27A.prop
date

