#
# Define settings for building locally
#
# Expects the environment variable NAGLIB_BASE_DIR to be set to the path of the
# base directory of local NAG library installation.
# 
# One of the environment variables MKLROOT and MKL_INCLUDE_DIR also needs to be
# set. MKLROOT if set should point to the root directory of the MKL library
# installation with `include`, `lib` etc. subdirectories. MKL_INCLUDE_DIR if
# set should point directly to the directory containing the MKL header files.

# Local library locations
ifdef MKL_INCLUDE_DIR
    LOCAL_MKL_DIR = $(MKL_INCLUDE_DIR)
else ifdef MKLROOT
    LOCAL_MKL_DIR=$(MKLROOT)/include
else
    $(error One of MKL_INCLUDE_DIR and MKLROOT needs to be set)
endif
ifdef NAGLIB_BASEDIR
    LOCAL_NAG_BLOCKS = $(NAGLIB_BASEDIR)/lp64/nag_interface_blocks
else
    $(error NAGLIB_BASEDIR needs to be set)
endif

# Local compiler
LOCAL_FC = gfortran

# Compiler and linker options
LOCAL_FL_LIBS = ${NAG_DIR}/lp64/lib/libnag_mkl.a -Wl,--no-as-needed \
                -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -liomp5 -lpthread \
                -lm -ldl

# Compiler flags
LOCAL_FC_FLAGS = -cpp -O3 -fcheck=bounds -ffree-line-length-0 -fimplicit-none \
                 -march=native -funroll-loops -m64 -fopenmp
LOCAL_FC_DEBUG_FLAGS = -cpp -g -fcheck=all -ffree-line-length-0 \
                       -fimplicit-none -fbacktrace -finit-real=nan \
                       -ffpe-trap=invalid,zero,overflow,underflow

# Local MPI settings
LOCAL_MPI_FLAGS = -DUSE_MPI=1
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC = mpifort