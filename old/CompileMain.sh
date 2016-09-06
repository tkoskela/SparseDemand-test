#!/bin/sh
#
#  script1 - Compile, link and execute low rank quadratic utility program
#
# Revision History

# 2015AUG12 LN  update to new intel, new NAG libraries
# 2015MAY24 :M  update path to INTL_MKL libraries
# 14aug2013 LN  update path to NAG libraries
# 07dec2012 LN   adapted from Horowitz
#
#
#  The shell variables must be set as follows:
#
#  NAG_DIR  :   the root directory of the NAG Fortran Library materials
#  fcompile :   the command (with options) to execute the Fortran compiler
#  flink    :   the linker option required to link the example program with
#               the NAG Fortran Library.
#
# Name of computer on which code is being compiled
HOST=`hostname`

# Clear the window
clear

# Choose which version to compile: DGB or EXE
#   script1 DBG    to compile debuggable
#   script1 EXE    to compile executable
#   script1        to compile default option: currently same as DBG
#   script1 dbg USE_MPI=1   to compile dbg and use MPI
#   script1 exe USE_MPI=1   to compile exe and use MPI
case $1 in
'DBG'|'dbg')
VERSION=$1
;;
'EXE'|'exe')
VERSION=$1
;;
*)
echo 'No valid version selected. Assuming default value: DBG.'
VERSION='DBG'
;;
esac


#  OPTIONS FOR UCL SERVER
case ${HOST} in (hpc* | node*) 
  echo $HOST
  NAG_DIR="/cm/shared/apps/NAG/fll6i25dc"  
  intel_mkl_dir="/cm/shared/apps/intel/composer_xe_2015.1.133/mkl/include"
  NAG_MKL_DIR="${NAG_DIR}/mkl_intel64_11.2.0/lib"
  NAG_LIBS="${NAG_DIR}/lib/libnag_mkl.a -Wl,--start-group ${NAG_MKL_DIR}/libmkl_core.a ${NAG_MKL_DIR}/libmkl_intel_lp64.a ${NAG_MKL_DIR}/libmkl_intel_thread.a -Wl,--end-group -liomp5 -lpthread -lm"

  fl_flags="${NAG_LIBS} -lifport -shared-intel"
  fc_exe_flags="-fpp -cxxlib -align all -funroll-loops -allow nofpp-comments -mtune=core2 -O3 -parallel -ipo -no-prec-div -xHost -axSSE4.2,SSSE3"
  fc_dbg_flags="-fpp -cxxlib -align all -allow nofpp_comments -g -debug extended -parallel";;
# options fo CS server
('vic.local')
  MPI_DIR=/share/apps/openmpi-1.8.1-intel-13.0.1/bin
  NAG_DIR=/share/apps/NAG/lib/fll6i24dcl
  intel_mkl_dir="/cm/shared/apps/intel/mkl/10.2.7.041/include"
  mkl_dir="${NAG_DIR}/mkl_intel64"
  nag_libs="${NAG_DIR}/lib/libnag_mkl.a -Wl,--start-group ${mkl_dir}/libmkl_core.a ${mkl_dir}/libmkl_intel_lp64.a ${mkl_dir}/libmkl_intel_thread.a -Wl,--end-group -liomp5 -lpthread"
  fl_flags="${nag_libs} -lifport -shared-intel"
  fc_exe_flags="-fpp -cxxlib -align all -funroll-loops -allow nofpp-comments -mtune=core2 -O3 -parallel -ipo -no-prec-div -xHost -axSSSE3"
  fc_dbg_flags="-fpp -cxxlib -align all -allow nofpp_comments -g -debug extended -parallel";;
#   OPTIONS FOR IFS SERVER
('penfold')
  NAG_DIR=/opt/NAG/fll6i21dcl      
  mkl_dir="${NAG_DIR}/mkl_em64t"   
  fl_flags="${NAG_DIR}/lib/libnag_mkl.a -L${mkl_dir} -lmkl_lapack -lmkl_em64t -lguide -lpthread -lirc"
  fc_exe_flags="-fpp -cxxlib-gcc -align all -allow nofpp_comments -g -debug extended"
  fc_dbg_flags="-fpp -ipo -O3 -xW -cxxlib-gcc -align all -parallel -unroll -allow nofpp_comments -diag-error-limit5 -DSLEEP=1";;
esac


case $2 in
'USE_MPI=1' | '1')
echo 'MPI version.'
fcompile="mpif90"
MPI_FLAGS="-debug parallel -DUSE_MPI=1"
;;
*)
echo 'Serial version (w/o MPI).'
fcompile="ifort"
MPI_FLAGS="-DUSE_MPI=0"
;;
esac

d1="include"
d2="source"
nagmod="-I${NAG_DIR}/nag_interface_blocks"
intelmkl="-I${intel_mkl_dir}"
mod1="${d1}/includeNR.f90 ${d1}/includeTools.f90 ${d1}/includePropertyList.f90"
mod2="${d2}/NewTools.f90 ${d2}/GlobalModule.f90 ${d2}/DataModule.f90 ${d2}/OutputModule.f90"
mod3="${d2}/LikelihoodModule.f90"
main="${d2}/main.f90"
out="SparseDemand"
case $2 in
'USE_MPI=1' | '1')
out="SparseDemand_mpi"
;;
esac

#
#clear

case ${VERSION} in ('DBG'|'dbg') 
echo "Debug version."
echo "Compiling, linking. In progress."
#   Use this line to compile, link and prepare for debugging.
${fcompile} ${fc_dbg_flags} ${MPI_FLAGS} ${nagmod} ${intelmkl} ${mod1} ${mod2} ${mod3} ${main} -o ${out}.dbg ${fl_flags}

;;
('EXE'|'exe')
echo "Optimized version."
echo "Compiling, linking. In progress."
#   Use this line to compile, link and prepare executable version.
${fcompile} ${fc_exe_flags} ${MPI_FLAGS} ${mod1} ${nagmod} ${intelmkl} ${mod2} ${mod3} ${main} -o ${out}.exe ${fl_flags}
;;
esac

echo "Compilation and linking completed."
#clear
