#
# Makefile - builds Numerical Recipes library
#
# To build a target, enter:
#
#   make <target>
#
# modification history
# --------------------
# 2021MAR06 LN edit to pathnames
# 13feb08 lpn cleanup local flags.
# 07feb08 bss adapted to be platform independent.
# 17aug07 bss added targets to handle recursive calling of Makefile.
# 06aug07 bss written.
#
#

.SUFFIXES: .a .o .f95 .f90 .mod

#---------------------------------------
# Define macros for Makefile

# Select Current Platform-specific parameters
WorkDir:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/../..
include $(WorkDir)/config/LoadLocalDefs.inc

# Compiler toolkit

# Fortran Compiler
FC = $(LOCAL_FC)

# Fortran Linker
FL = $(LOCAL_FL)

# Archive tools -- use GNU by default

AR = $(LOCAL_AR)
RANLIB = $(LOCAL_RANLIB)
AR_FLAGS = $(LOCAL_AR_FLAGS)

# Where libraries live
NR_DIR = $(LOCAL_NR_DIR)
MKL_DIR = ${LOCAL_MKL_DIR}

# Compiler and linker options
FL_LIBS = $(LOCAL_FL_LIBS) -I${MKL_DIR}

# Fortran compiler flags
#
# Note: -fast == -ipo -O3 -no-prec-div -xP -static
#       -static causes link errors
#       -no-prec-div reduces accuracy and increases speed
# Note: -cxxlib-gcc is needed by NAG to prevent warning about glibc
#

# Note: optimization removes symbols from the .o.
FC_FLAGS = $(LOCAL_FC_FLAGS)
FC_DEBUG_FLAGS = $(LOCAL_FC_DEBUG_FLAGS)

#---------------------------------------
# Source code for numerical recipies library
NR_SRC = \
	nrtype.f90 \
	nrutil.f90 \
	nr.f90 \
	ran_state.f90 \
	numerical_recipes.f90

OBJS = $(patsubst %.f90, %.o, $(NR_SRC) )
DEBUG_OBJS = $(patsubst %.f90, %_d.o, $(NR_SRC) )

#---------------------------------------
# Targets - these are objects which can be built
#

all: clean libnr.a libnr_d.a

exe: libnr.a

dbg: libnr_d.a

# Rule to build libraries

# The stuff with $(patsubst ...) below is 'syntactic sugar' to
# prefix each source file with the path to where it is stored
libnr.a : $(OBJS)
#	$(FC) -fpic -c $^
	$(AR) $(AR_FLAGS) $@ $(OBJS)
	$(RANLIB) $@

libnr_d.a : $(DEBUG_OBJS)
#	$(FC) -fpic $(FC_DEBUG_FLAGS) -c $^
#	$(FC) -fpic -c $^
	$(AR) $(AR_FLAGS) $@ $(DEBUG_OBJS)
	$(RANLIB) $@

# clean up any messes which have been created
clean:
	rm -f *.o *.mod a.out *.a *.dbg *.exe


# Pattern Rules

$(OBJS) : %.o : %.f90
	$(FC) -fpic $(FC_FLAGS) -c $< -o $@

$(DEBUG_OBJS) : %_d.o : %.f90
	$(FC) -fpic $(FC_DEBUG_FLAGS) -c $< -o $@
