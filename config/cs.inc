#
# config.cs.inc - defines local settings for CS cluster
#
# modification history
# --------------------
# 20191112  LN  edit for new cs cluster
# 2018JUL31 LN  updated path to intel 2018 and NAG Mark 26
# 2017JAN03 LN  clean up config file
# 13aug2015 LN  update path to NAG
# 14aug2013 LN  update path to NAG
# 19nov2012 LN  update NAG flags and Intel flags
#
#-------------------
# Local MPI settings
# mpif90 sets all other flags.
LOCAL_MPI_DIR         = /share/apps/intel_parallel_studio_xe_2018/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin
LOCAL_MPI_LIBS        =
LOCAL_MPI_FLAGS       = -DUSE_MPI=1 
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC          = mpiifort
LOCAL_MPI_CC          = mpicc

#-------------------
# Local NAG settings

# Library locations
LOCAL_MKL_DIR     = /share/apps/intel_parallel_studio_xe_2018/compilers_and_libraries_2018.3.222/linux/mkl/include

# fll6i26dcl:   intel fortran compiler version
# fll6i26dfl:   gfortran compiler version
LOCAL_NAG_DIR     = /share/apps/NAG/Fortran/fll6i26dcl
LOCAL_NAG_BLOCKS  = ${LOCAL_NAG_DIR}/nag_interface_blocks
LOCAL_KNITRO_DIR  = /share/apps/econ/knitro/10.0.1-z
LOCAL_KNITRO_LIBS = $(LOCAL_KNITRO_DIR)/lib/libknitro.so
KNITRORELEASE = 1001
USE_NAG = 1

#-------------------
# Local Compiler Toolkit settings

LOCAL_CC       = icc
LOCAL_FC       = ifort
LOCAL_FL       = ifort
LOCAL_AR       = xiar
LOCAL_AR_FLAGS = rc
LOCAL_RANLIB   = ranlib

# Compiler and linker options
LOCAL_FL_LIBS = $(LOCAL_NAG_DIR)/lib/libnag_mkl.a \
                -Wl,--start-group \
                ${LOCAL_NAG_DIR}/mkl_intel64_11.3.3/lib/libmkl_core.a \
                ${LOCAL_NAG_DIR}/mkl_intel64_11.3.3/lib/libmkl_intel_lp64.a \
                ${LOCAL_NAG_DIR}/mkl_intel64_11.3.3/lib/libmkl_intel_thread.a -Wl,--end-group \
                -liomp5 -lpthread -lm -shared-intel \
                -ldl -lstdc++
LOCAL_CC_LIBS = -ldl -lstdc++

# Fortran compiler flags
#
# Note: -fast == -ipo -O3 -no-prec-div -xP -static
#        -static causes link errors
#        -no-prec-div leads to lower accuracy and more speed.
# Note: -cxxlib is needed by NAG to prevent warning about glibc

# CS Compiler options
LOCAL_CC_FLAGS = -ipo -O3
# -no-prec-div       faster but not full precision for division
# -fp-model fast=2   faster but not full precision
# -mtune=???
# -march=
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -assume realloc_lhs \
                 -prec-div -fp-model precise -ipo \
                 -funroll-loops -parallel -O3 -xAVX 
#LOCAL_FC_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -assume realloc_lhs \
#                 -prec-div -fp-model precise -ipo \
#                 -funroll-loops -parallel -O3 -xHost 

LOCAL_FC_DEBUG_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -assume realloc_lhs \
                       -g -debug extended -gdwarf-2 -check bounds
LOCAL_CC_DEBUG_FLAGS = -fopenmp -cxxlib -align -g -debug extended 

