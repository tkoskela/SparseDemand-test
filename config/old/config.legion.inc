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
# 23JAN2014 LN  edit to work on legion, adapt from ucl.config.inc
# 14aug2013 LN  update path to NAG
# 14aug2013 LN  update path to NAG
# 19nov2012 LN  update NAG flags and Intel flags
# 16jun2011 LPN update compiler flags
# 01sep2009 LPN New flags to link to NAG libraries. Old flags don't work with new license.
# 13feb2008 lnp edited mpi flags
# 08feb2008 bss added flag to fix warning about feupdateenv.
# 07feb2008 bss added MPI support.
# 06feb2008 bss written.
#
#-------------------
# Local MPI settings
# mpif90 sets all other flags.
LOCAL_MPI_LIBS         =
LOCAL_MPI_FLAGS       = -DUSE_MPI=1 
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC          = mpif90

#-------------------
# Local NAG settings

# Where libraries live
LOCAL_NAG_DIR = /shared/ucl/apps/NAG/fortran-library/fll6i24dcl
LOCAL_NR_DIR = /home/uctpln0/GeneralCode/nr
LOCAL_TOOLS_DIR = /home/uctpln0/GeneralCode/tools

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
                ${LOCAL_NAG_DIR}/mkl_intel64/libmkl_core.a \
                ${LOCAL_NAG_DIR}/mkl_intel64/libmkl_intel_lp64.a \
                ${LOCAL_NAG_DIR}/mkl_intel64/libmkl_intel_thread.a -Wl,--end-group \
                -liomp5 -lpthread                        \
                -lifport                         \
                -shared-intel

# Fortran compiler flags
#
# Note: -fast == -ipo -O3 -no-prec-div -xP -static
#        -static causes link errors
#        -no-prec-div leads to lower accuracy and more speed.
# Note: -cxxlib is needed by NAG to prevent warning about glibc
# to improve performance
# -ax<processor>
# -ftz  flush denormal to zero (not needed if -O3 is set)
# -march=<processor>   (core-avx2 advanced vector extensions 2)
#
# -mtune=<processor>
# -O[n]
# -prof-gen
# -prof-use
# -x<processor>
# fast = -ipo -O3 -no-prec-div  -static -xHost
# allowed processors are:
#  test-node = SSE4.2 POPCNT instructions
# questions
#  -ax<???>
#  -x???   (CORE-AVX2, CORE-AVX-I, AVX, SSE4.2, SSE4.1, SSSE3,SSE3,SSE2,Host)
#  -march=<???>
#  -mtune=<???>
# UCL_HPC Compiler options
#LOCAL_FC_FLAGS = -fpp -cxxlib -align all -funroll-loops \
#                 -allow nofpp_comments \
#                 -mtune=core2 -O3 -parallel -ipo -no-prec-div -xSSE4.2 -axSSE4.2
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -funroll-loops \
                 -allow nofpp_comments \
                 -O3 -parallel -ipo -no-prec-div -mtune=core2 -xSSSE3
LOCAL_FC_DEBUG_FLAGS = -fpp -cxxlib -align all -allow nofpp_comments -g -debug extended -parallel

