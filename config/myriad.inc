#
# config.myriad.inc - defines local settings for myriad cluster
#
# modification history
# --------------------
# 2021MAR06 LN  edit path names
# 2019DEC15 LN  edit for myriad
# 2019AUG29 LN  edit for myriad cluster
# 23JAN2014 LN  edit to work on legion, adapt from ucl.config.inc
#
# On myriad use modules to set paths to libraries.
# tested with the following module environment (TK March 2021):
# module purge
# module load default-modules/2018
# module load compilers/intel/2018/update3
# module load nag/fortran/mark26/intel-2017
# module load knitro/12.0.0/gnu-4.9.2
#
#-------------------
# Local MPI settings
# mpif90 sets all other flags.
LOCAL_MPI_DIR         = $(MPI_HOME)/intel64/bin
LOCAL_MPI_LIBS        =
LOCAL_MPI_FLAGS       = -DUSE_MPI=1
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC          = mpiifort
LOCAL_MPI_CC          = mpicc

#-------------------
# Local NAG settings

# Library locations
LOCAL_MKL_DIR     = $(MKLROOT)/include
LOCAL_NAG_BLOCKS  = $(NAG_DIR)/nag_interface_blocks
LOCAL_KNITRO_DIR  = $(KNITRODIR)
LOCAL_KNITRO_LIBS = $(LOCAL_KNITRO_DIR)/lib/libknitro.so
KNITRORELEASE = 1200
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
LOCAL_FL_LIBS = -lnag_mkl -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -shared-intel -ldl -lstdc++
LOCAL_CC_LIBS = -ldl -lstdc++

# Fortran compiler flags
#
# Note: -fast == -ipo -O3 -no-prec-div -xP -static
#        -static causes link errors
#        -no-prec-div leads to lower accuracy and more speed.
# Note: -cxxlib is needed by NAG to prevent warning about glibc

# MYRIAD Compiler options
LOCAL_CC_FLAGS = -ipo -O3
# -no-prec-div       faster but not full precision for division
# -fp-model fast=2   faster but not full precision
# -mtune=core2
# -xcode=icelake-server -mtune=icelake-server -march=cascadelake
# -mtune=core-avx2
# -march=core-avx2
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -assume realloc_lhs \
                 -prec-div -fp-model precise -ipo -funroll-loops -parallel -O3 \
                 -xHost

LOCAL_FC_DEBUG_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -assume realloc_lhs \
                       -g -debug extended -gdwarf-2 -check bounds
LOCAL_CC_DEBUG_FLAGS = -fopenmp -cxxlib -align -g -debug extended
