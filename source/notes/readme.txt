Test version of SparseDemand

# Preparation
1) Load required modules
   module load nag/fortran/mark26/intel-2017
2) Create output directory.
3) Edit name of OutDir in ../inputs/Test1/A1a.estimate.prop 
   to match directory created in 2). 

# Compile and run serial version 
1) compile exe
   make cleanall
   make exe
2) Run
   ./SparseDemand.exe ../inputs/Test1/A1a.estimate.prop

  a) std. io. shows iterations of optimization software.
  b) simulated data and estimation results are saved to OutDir.

# Compile and run parallel version
1) compile exe
   rm *0.0 *.mod
   make USE_MPI=1 exe
2) Run
   mpirun -np 8 ./SparseDemand_mpi.exe ../inputs/Test1/A1a.estimate.prop

  a) std. io. shows iterations of optimization software.
  b) simulated data and estimation results are saved to OutDir.
  c) note results are not identical to serial version because the random numbers are generated from different seed


 
