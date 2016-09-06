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
# 12aug2015 LN updated NAG path to new version
# 14aug2013 LN updated NAG path to new version
# 19nov2012 LN  update NAG and Intel flags
# 30jun2011 LN  updated flags to new NAG libraries and new Intel compiler
# 01sep2009 LPN New flags to link to NAG libraries. Old flags no longer work
# 13feb2008 lpn cleanup local flags, etc.
# 07feb2008 bss added MPI support.
# 06feb2008 bss written.
#
#-------------------
# Local MPI settings
# No special flags required for mpi. all are set by mpif90.

LOCAL_MPI_LIBS        =
LOCAL_MPI_FLAGS       = -DUSE_MPI=1
LOCAL_MPI_DEBUG_FLAGS = -DUSE_MPI=1
LOCAL_MPI_FC          = mpif90

#-------------------
# Local NAG settings

# Where libraries live
LOCAL_NAGDIR = /cm/shared/apps/NAG/fll6i25dc
LOCAL_NAG_MKL_DIR = $(LOCAL_NAGDIR)/mkl_intel64_11.2.0/lib
LOCAL_NR_DIR = ../nr
LOCAL_TOOLS_DIR = ../tools
LOCAL_MKL_DIR = /cm/shared/apps/intel/composer_xe_2015.1.133/mkl/include

#-------------------
# Local Complier Toolkit settings

LOCAL_FC       = ifort
LOCAL_FL       = ifort
LOCAL_AR       = xiar
LOCAL_AR_FLAGS = rc
LOCAL_RANLIB   = ranlib

# Compiler and linker options
LOCAL_FL_LIBS = $(LOCAL_NAG_DIR)/lib/libnag_mkl.a \
                -Wl,--start-group \
                ${LOCAL_NAG_MKL_DIR}/libmkl_core.a \
                ${LOCAL_NAG_MKL_DIR}/libmkl_intel_lp64.a \
                ${LOCAL_NAG_MKL_DIR}/libmkl_intel_thread.a -Wl,--end-group \
                -liomp5 -lpthread -lm                 \
                -lifport                         \
                -shared-intel

# Fortran compiler flags
#
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -funroll-loops \
                 -allow nofpp_comments \
                 -mtune=core2 -O3 -parallel -ipo -no-prec-div -xHost -axSSSE3
LOCAL_FC_DEBUG_FLAGS =	-fpp -g -debug extended -cxxlib \
                        -allow nofpp_comments -parallel

