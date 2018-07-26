#
# Makefile - builds Sparse demand model
#
# To build a target, enter:
#
#   make <target>
#
# Targets:
#
#   all	  	 - builds everything
#   clean	 - deletes all .o, .a, binaries, etc.
#   nuke     - clean + clean all dependencies
#   exe		 - optimized application
#   dbg		 - application with debug symbols
#
# Options:
#
#   USE_MPI	- defaults to 0 (no MPI).  If defined to 1, build 
#                 application with MPI libraries for execution on 
#                 multiple processors
#
# Example:
#
#   To build an optimized application with MPI:
#
#     make USE_MPI=1 exe
#
#   To build a optimized application for one processor (sequential):
#
#     make exe
#
#   To build a sequential debug application:
#
#     make dbg
#
# modification history
# --------------------
# 2018JUL25 LN   clean up and edit
# 09OCT2013 LN   load ApproximationModule.o before ValueModule.o since the latter uses the former
# 29jun2011 LPN  rename *_module to *Module
# 11nov2010 LPN  adapt from housing project

.SUFFIXES: .a .o .f95 .f90 .mod

# Select correct Macros based on platform in use
include ../config/LoadLocalDefs.inc

# Compiler toolkit

# Fortran Compiler
FC = $(LOCAL_FC)

# Default is not to use MPI libraries
# USE_MPI = 1 --> build for use with MPI

ifeq ($(USE_MPI),1)
  LOCAL_FL_LIBS += $(LOCAL_MPI_LIBS)
  LOCAL_FC_FLAGS += $(LOCAL_MPI_FLAGS)
  LOCAL_FC_DEBUG_FLAGS += $(LOCAL_MPI_DEBUG_FLAGS) -DSLEEP=1
  FC = $(LOCAL_MPI_FC)
endif

# Archive tools -- use GNU by default

AR = $(LOCAL_AR)
RANLIB = $(LOCAL_RANLIB)
AR_FLAGS = $(LOCAL_AR_FLAGS)

# Where libraries live
NR_DIR = $(LOCAL_NR_DIR)
TOOLS_DIR = $(LOCAL_TOOLS_DIR)
MKL_DIR = ${LOCAL_MKL_DIR}
NAG_DIR = ${LOCAL_NAG_DIR}
OUT_DIR = ../lib

# Compiler and linker options
FL_LIBS = $(LOCAL_FL_LIBS) 
FC_FLAGS = $(LOCAL_FC_FLAGS) -I${MKL_DIR} -I$(NR_DIR) -I$(TOOLS_DIR) -I${NAG_DIR}/nag_interface_blocks -DHOST_ID=$(HOST_ID)

# Add the following line to FC_FLAGS to use profiling data to 
# improve compiler optimization:
#		 -prof-use -prof-dir ./profdata
# 
# Add the following line to FC_FLAGS to generate profile data
#
#		 -prof-genx -prof-dir ./profdata

FC_DEBUG_FLAGS = $(LOCAL_FC_DEBUG_FLAGS) -I${MKL_DIR} -I$(NR_DIR) -I$(TOOLS_DIR) -I${NAG_DIR}/nag_interface_blocks -DHOST_ID=$(HOST_ID) 

#---------------------------------------
# Source code for numerical recipes library
NR_SRC = \
	nrtype.f90 \
	nrutil.f90 \
	nr.f90 \
	ran_state.f90 \
	numerical_recipes.f90

TOOLS_SRC = \
	ToolsModule.f90 \
	PropertyList.f90

# Modules -- not the main program!
MOD_OBJS = \
	NewTools.o \
	GlobalModule.o \
	OutputModule.o \
	DataModule.o \
	LikelihoodModule.o 

# The main program
PRG_OBJS = main.o

OBJS = $(MOD_OBJS) $(PRG_OBJS)
DEBUG_OBJS = $(patsubst %.o, %_d.o, $(OBJS))

# Application Name -- what to build
ifeq ($(USE_MPI),1)
TGT = SparseDemand_mpi
else
TGT = SparseDemand
endif

# NOTE:  This order of single and double quotes is crucial
# because any other order will cause mpiifort to fail because
# its preprocessor is weaker than ifort's.  In particular, it
# seems to only take one pass and have trouble with quoting...
FC_FLAGS       += -DAPP_NAME='"$(TGT)"'
FC_DEBUG_FLAGS += -DAPP_NAME='"$(TGT)"'

#---------------------------------------
# Targets - these are objects which can be built
#

exe:    $(TGT).exe

dbg:    $(TGT).dbg

all: clean dbg exe

nuke:	cleantools cleannr clean

#-----------------------------
# Rule to build libraries

# The codes $(patsubst ...) adds a prefix to 
# each source file defining the path where it is stored
# This section also builds the libraries (libnr,libtools) if they have not yet already been built.
#-----------------------------

$(NR_DIR)/libnr.a : $(patsubst %, $(NR_DIR)/%, $(NR_SRC) )
	cd $(NR_DIR) ; make libnr.a

$(NR_DIR)/libnr_d.a : $(patsubst %, $(NR_DIR)/%, $(NR_SRC) )
	cd $(NR_DIR) ; make libnr_d.a

$(TOOLS_DIR)/libtools.a : $(patsubst %, $(TOOLS_DIR)/%, $(TOOLS_SRC) )
	cd $(TOOLS_DIR) ; make libtools.a

$(TOOLS_DIR)/libtools_d.a : $(patsubst %, $(TOOLS_DIR)/%, $(TOOLS_SRC) )
	cd $(TOOLS_DIR) ; make libtools_d.a

#-----------------------------
# Rules to build application
$(TGT).exe : $(NR_DIR)/libnr.a $(TOOLS_DIR)/libtools.a $(OBJS)
	$(FC) $(FC_FLAGS) $^ -o $@	\
		$(FL_LIBS)		\
		$(NR_DIR)/libnr.a $(TOOLS_DIR)/libtools.a

$(TGT).dbg : $(NR_DIR)/libnr_d.a $(TOOLS_DIR)/libtools_d.a $(DEBUG_OBJS)
	$(FC) $(FC_DEBUG_FLAGS) $^ -o $@ \
		$(NR_DIR)/libnr_d.a $(TOOLS_DIR)/libtools_d.a	\
		$(FL_LIBS)

#-----------------------------
# clean up any messes which have been created
clean:
	rm -f *.o $(TGT) *.mod a.out *.a *.dbg *.exe

cleantools:
	cd $(TOOLS_DIR) ; make clean

cleannr:
	cd $(NR_DIR) ; make clean

# Pattern Rules

$(OBJS) : %.o : %.f90 
	$(FC) $(FC_FLAGS) -c $< -o $@

$(DEBUG_OBJS) : %_d.o : %.f90 
	$(FC) $(FC_DEBUG_FLAGS) -c $< -o $@

