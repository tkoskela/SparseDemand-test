# make dbg                compile dbg version of library
# make TestinverseNormal  compile test of inverse normal
#
# modification history
# --------------------
# 2021MAR06 LN edit path names 
# 2018AUG20 LN update to work with Intel 2018
#              add option to compile Test InverseNormal

.SUFFIXES: .o .f95 .f90 .mod .a

# Select Corrent Platform-specific parameters
include ../../config/LoadLocalDefs.inc

# Fortran Compiler
FC = $(LOCAL_FC)

# Fortran Linker
FL = $(LOCAL_FL)

# Archive tools -- use GNU by default
AR = $(LOCAL_AR)
RANLIB = $(LOCAL_RANLIB)
AR_FLAGS = $(LOCAL_AR_FLAGS)

# Compiler and linker options
FL_LIBS = $(LOCAL_FL_LIBS)

# Numerical Recipes library
NR_DIR = $(LOCAL_NR_DIR)
MKL_DIR = ${LOCAL_MKL_DIR}
NAG_INTERFACE_BLOCKS = $(LOCAL_NAG_BLOCKS)

# Fortran compiler flags
#

# Note: optimization eliminates symbols from the .o...
FC_FLAGS = $(LOCAL_FC_FLAGS) -I$(NR_DIR) -I${MKL_DIR} -I$(NAG_INTERFACE_BLOCKS)
FC_DEBUG_FLAGS = $(LOCAL_FC_DEBUG_FLAGS) -I$(NR_DIR) -I${MKL_DIR} -I$(NAG_INTERFACE_BLOCKS)

#---------------------------------------
# Source for Tools library

TOOLS_SRC = \
	ToolsModule.f90 \
	SplineTools.f90 \
	MultiModule.f90 \
	PropertyList.f90

OBJS = $(patsubst %.f90, %.o, $(TOOLS_SRC) )
DEBUG_OBJS = $(patsubst %.f90, %_d.o, $(TOOLS_SRC) )

#---------------------------------------
# Targets - these are objects which can be built
#

all: clean libtools.a libtools_d.a

exe: libtools.a

dbg: libtools_d.a

# Rule to build libraries

# The stuff with $(patsubst ...) below is 'syntactic sugar' to
# prefix each source file with the path to where it is stored

libtools.a : $(NR_DIR)/libnr.a $(OBJS)
	$(AR) $(AR_FLAGS) $@ $(OBJS)
	$(RANLIB) $@

libtools_d.a : $(NR_DIR)/libnr_d.a $(DEBUG_OBJS)
	$(AR) $(AR_FLAGS) $@ $(DEBUG_OBJS)
	$(RANLIB) $@

$(NR_DIR)/libnr.a:
	make -C $(NR_DIR) libnr.a

$(NR_DIR)/libnr_d.a:
	make -C $(NR_DIR) libnr_d.a

# clean up any messes which have been created
clean:
	rm -f *.o *.mod a.out *.a *.dbg *.exe

# Pattern Rules

$(OBJS) : %.o : %.f90
	$(FC) -fpic $(FC_FLAGS) -c $< -o $@

$(DEBUG_OBJS) : %_d.o : %.f90
	$(FC) -fpic $(FC_DEBUG_FLAGS) -c $< -o $@
