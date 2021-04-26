#
# Define settings for building locally
#

# Local MPI settings
LOCAL_MPI_DIR = $(MPI_HOME)/intel64/bin
LOCAL_MPI_FLAGS = -DUSE_MPI=1
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC = mpifort
LOCAL_MPI_CC = mpicc

# Local library locations
LOCAL_MKL_DIR = $(MKL_INCLUDE_DIR)
LOCAL_NAG_BLOCKS = $(NAGLIB_BASEDIR)/lp64/nag_interface_blocks
LOCAL_NR_DIR = $(WorkDir)/GeneralCode/nr
LOCAL_TOOLS_DIR = $(WorkDir)/GeneralCode/tools
LOCAL_KNITRO_DIR = $(KNITRODIR)
LOCAL_KNITRO_LIBS = $(LOCAL_KNITRO_DIR)/lib/libknitro.so
KNITRORELEASE = 1200
USE_NAG = 1

# Local Compiler Toolkit settings
LOCAL_CC = gcc
LOCAL_FC = gfortran
LOCAL_FL = gfortran
LOCAL_AR = ar
LOCAL_AR_FLAGS = rc
LOCAL_RANLIB = ranlib

# Compiler and linker options
LOCAL_FL_LIBS = ${NAG_DIR}/lp64/lib/libnag_mkl.a  -Wl,--no-as-needed \
                -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -liomp5 -lpthread -lm -ldl
LOCAL_CC_LIBS = -ldl -lstdc++


# Fortran compiler flags
LOCAL_FC_FLAGS = -cpp -fcheck=bounds -O3  -march=native -funroll-loops -m64 -fopenmp
LOCAL_FC_DEBUG_FLAGS = -cpp -g -Wall -Wextra -Wconversion -fcheck=bounds -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan

# C compiler flags
LOCAL_CC_FLAGS = -ipo -O3
LOCAL_CC_DEBUG_FLAGS = -fopenmp -cxxlib -align -g -debug extended
