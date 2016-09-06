#
# config.ifs.hpc.inc - defines local settings for IFS HPC environment
#
# You should define the local version of any platform-specific
# settings or configuration elements here.  E.g., the location
# of the NAG libraries or MPI.  You should not put platform-specific
# information in the Makefile!!!
#
# modification history
# --------------------
# 13DEC2014 LN  update NAG library path
# 13DEC2014 LN  update NAG library path
# 13DEC2014 LN  update NAG library path
# 13MAR2013 LN updated to IFS HPC
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
LOCAL_MPI_DEBUG_FLAGS = -debug parallel -DUSE_MPI=1
LOCAL_MPI_FC          = mpiifort

#-------------------
# Local NAG settings

# Where libraries live
LOCAL_NAGDIR = /opt/NAG/fll6i24dcl
LOCAL_NR_DIR = ../nr
LOCAL_TOOLS_DIR = ../tools

#-------------------
# Local Compiler Toolkit settings

LOCAL_FC       = ifort
LOCAL_FL       = ifort
LOCAL_AR       = xiar
LOCAL_AR_FLAGS = rc
LOCAL_RANLIB   = ranlib

# Compiler and linker options
# Flags to link to NAG
LOCAL_FL_LIBS =	$(LOCAL_NAGDIR)/lib/libnag_mkl.a \
		            -Wl,--start-group \
		            ${LOCAL_NAGDIR}/mkl_intel64/libmkl_core.a \
                ${LOCAL_NAGDIR}/mkl_intel64/libmkl_intel_lp64.a \
                ${LOCAL_NAGDIR}/mkl_intel64/libmkl_intel_thread.a -Wl,--end-group \
		            -liomp5 -lpthread \
		            -lifport \
                -shared-intel
# OLD flags: replaced 2012NOV19
#LOCAL_FL_LIBS =	$(LOCAL_NAGDIR)/lib/libnag_mkl.a \
#		            -L$(LOCAL_NAGDIR)/mkl_em64t -Wl,--start-group \
#		            ${LOCAL_NAGDIR}/mkl_em64t/libmkl_intel_lp64.a \
#                ${LOCAL_NAGDIR}/mkl_em64t/libmkl_intel_thread.a \
#                ${LOCAL_NAGDIR}/mkl_em64t/libmkl_core.a -Wl,--end-group \
#		            -lpthread \
#		            -lifport \
#                -shared-intel

# Fortran compiler flags
#
LOCAL_FC_FLAGS = -fpp -cxxlib -align all -funroll-loops \
                 -allow nofpp_comments \
                 -mtune=core2 -O3 -parallel -ipo -no-prec-div -xHost -axSSSE3
LOCAL_FC_DEBUG_FLAGS =	-fpp -cxxlib -align all -allow nofpp_comments -g -debug extended -parallel

