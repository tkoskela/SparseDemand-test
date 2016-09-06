#
# config.cs.inc - defines local settings for CS cluster
#
# You should define the local version of any platform-specific
# settings or configuration elements here.  E.g., the location
# of the NAG libraries or MPI.  You should not put platform-specific
# information in the Makefile!!!
#
# modification history
# --------------------
# 2015SEP15 LN  update paths for Intel and for NAG Mark 25
# 07OCT2014 LPN adapted from config.ucl.hpc.inc
# 14aug2013 LN  update path to new NAG libraries
# 24feb2013 LPN update to new NAG libraries
# 20jan2011 LPN update path to NAG 
# 16sep2009 LPN update NAG linking flags.
# 13feb2008 lnp edited mpi flags
# 08feb2008 bss added flag to fix warning about feupdateenv.
# 07feb2008 bss added MPI support.
# 06feb2008 bss written.
#
#-------------------
# Local MPI settings
# mpif90 sets all other flags.
LOCAL_MPI_DIR         = /share/apps/openmpi-1.10.1-intel-15.0.1/bin
LOCAL_MPI_LIBS        =
LOCAL_MPI_FLAGS       = -DUSE_MPI=1 
LOCAL_MPI_DEBUG_FLAGS = -DUSE_MPI=1
LOCAL_MPI_FC          = $(LOCAL_MPI_DIR)/mpifort

#-------------------
# Local NAG settings

# Where libraries live
LOCAL_NAG_DIR = /share/apps/NAG/lib/fll6i25dcl
NAG_MKL_DIR   = ${LOCAL_NAG_DIR}/mkl_intel64_11.2.0/lib
LOCAL_NR_DIR  = /home/uctpln0/GeneralCode/nr
LOCAL_TOOLS_DIR = /home/uctpln0/GeneralCode/tools

#-------------------
# Local Compiler Toolkit settings

LOCAL_MKL_DIR   = /share/apps/intel/composer_xe_2015.1.133/mkl/include
LOCAL_INTEL_DIR = /share/apps/intel/composer_xe_2015.1.133/bin/intel64
LOCAL_FC        = $(LOCAL_INTEL_DIR)/ifort
LOCAL_FL        = $(LOCAL_INTEL_DIR)/ifort
LOCAL_AR        = $(LOCAL_INTEL_DIR)/xiar
LOCAL_AR_FLAGS  = rc
LOCAL_RANLIB    = ranlib

# Compiler and linker options
LOCAL_FL_LIBS = $(LOCAL_NAG_DIR)/lib/libnag_mkl.a \
                -Wl,--start-group \
                $(NAG_MKL_DIR}/libmkl_core.a \
                $(NAG_MKL_DIR)/libmkl_intel_lp64.a \
                $(NAG_MKL_DIR)/libmkl_intel_thread.a -Wl,--end-group \
                -liomp5 -lpthread -lm -lifport -shared-intel

# Fortran compiler flags
#
# Note: -cxxlib is needed by NAG to prevent warning about glibc

# CS_HPC Compiler options
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -funroll-loops \
                 -allow nofpp_comments \
                 -O3 -parallel -ipo -no-prec-div \
                 -xHost -mtune=core2
LOCAL_FC_DEBUG_FLAGS = -fpp -g -O0 -debug extended -cxxlib \
		       -allow nofpp_comments -parallel

