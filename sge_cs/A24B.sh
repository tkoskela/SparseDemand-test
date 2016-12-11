#$ -S /bin/bash
#$ -wd /home/uctpln0/FuirtDemand/code/fortran/source
#$ -o /home/uctpln0/FruitDemand/code/fortran/output/A24B_mpi.out
#$ -N A24B
#$ -j y
#$ -cwd
#$ -l tmem=4G,h_vmem=4G
#$ -l h_rt=24:0:0
#$ -l hostname=burns*
#$ -R y
#$ -pe orte 128

cd /home/uctpln0/FruitDemand/code/fortran/source
source /home/uctpln0/GeneralCode/config/SetEnv.sh
. /share/apps/econ/Modules/modules.sh
module load gcc/5.2.0 nag/fll6i25dc intel/composer/2015.1.133 openmpi/intel/1.10.0 totalview knitro
date
#./SparseDemand.dbg ../inputs/A1.prop
#mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A31_5.prop
#mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A24/A24A.prop
#mpirun --mca btl_tcp_if_include 128.41.96.0/21 -np $NSLOTS SparseDemand_mpi.exe ../inputs/A24/A24B.prop
which mpirun
mpirun --mca btl_tcp_if_include 10.0.0.0/8 --mca mpi_preconnect_mpi 1 -np $NSLOTS SparseDemand_mpi.exe ../inputs/A24/A24B.prop
#mpirun --mca btl_tcp_if_include eth1 -np $NSLOTS SparseDemand_mpi.exe ../inputs/A24/A24B.prop
date

# --mca btl_base_verbose 30
# btl_if_exclude = lo
#--mca orte_base_help_aggregate 0 
# btl tcp,self
