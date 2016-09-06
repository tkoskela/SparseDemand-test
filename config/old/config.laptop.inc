#
# config.ucl.hpc.inc - defines local settings for UCL HPC environment
#
# You should define the local version of any platform-specific
# settings or configuration elements here.  E.g., the location
# of the NAG libraries or MPI.  You should not put platform-specific
# information in the Makefile!!!
#
# modification history
# --------------------
# 13DEC2014 LPN edit for laptop configuration
#
#-------------------
# Local MPI settings
# mpif90 sets all other flags.
LOCAL_MPI_LIBS        =
LOCAL_MPI_FLAGS       = -DUSE_MPI=0 
LOCAL_MPI_DEBUG_FLAGS = -DUSE_MPI=0
#LOCAL_MPI_FC          = mpiifort
LOCAL_MPI_FC          = gfortran

#-------------------
# Local NAG settings

# Where libraries live
LOCAL_NAGDIR = /cm/shared/apps/NAG/fll6i24dc

#-------------------
# Local Compiler Toolkit settings

LOCAL_FC       = gfortran
LOCAL_FL       = gfortran
LOCAL_AR       = xiar
LOCAL_AR_FLAGS = rc
LOCAL_RANLIB   = ranlib

# Compiler and linker options
LOCAL_FL_LIBS = ${LOCAL_NAGDIR}/lib/libnag_mkl.a \
                -Wl,--start-group \
                ${LOCAL_NAGDIR}/mkl_intel64/libmkl_core.a \
                ${LOCAL_NAGDIR}/mkl_intel64/libmkl_intel_lp64.a \
                ${LOCAL_NAGDIR}/mkl_intel64/libmkl_intel_thread.a -Wl,--end-group \
                -liomp5 -lpthread -lifport -shared-intel

# Fortran compiler flags
#
# Note: -fast == -ipo -O3 -no-prec-div -xP -static
#        -static causes link errors
#        -no-prec-div leads to lower accuracy and more speed.
# Note: -cxxlib is needed by NAG to prevent warning about glibc

# UCL_HPC Compiler options
LOCAL_FC_FLAGS = -cpp 
LOCAL_FC_DEBUG_FLAGS = -cpp -g -debug extended

