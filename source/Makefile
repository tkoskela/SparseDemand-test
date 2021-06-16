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
#   exe		 - optimized application
#   dbg		 - application with debug symbols
#   tests    - application to run tests
#
# Options:
#
#   USE_MPI	- defaults to 0 (no MPI).  If defined to 1, build application with
#             MPI libraries for execution on multiple processors
#
# Examples:
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

# Select correct variable definitions based on platform in use
# config directory assumed to be in parent of directory the makefile is in
WORK_DIR = $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/..
include $(WORK_DIR)/config/LoadLocalDefs.inc

# Default is not to use MPI libraries
# USE_MPI = 1 --> build for use with MPI
ifeq ($(USE_MPI),1)
  TGT = SparseDemand_mpi
  FC = $(LOCAL_MPI_FC)
  FL_LIBS = $(LOCAL_FL_LIBS) $(LOCAL_MPI_LIBS)
  FC_FLAGS = $(LOCAL_FC_FLAGS) $(LOCAL_MPI_FLAGS)
  FC_DEBUG_FLAGS = $(LOCAL_FC_DEBUG_FLAGS) $(LOCAL_MPI_DEBUG_FLAGS)
else
  TGT = SparseDemand
  FC = $(LOCAL_FC)
  FL_LIBS = $(LOCAL_FL_LIBS)
  FC_FLAGS = $(LOCAL_FC_FLAGS)
  FC_DEBUG_FLAGS = $(LOCAL_FC_DEBUG_FLAGS)
endif

FC_FLAGS += -I${LOCAL_MKL_DIR} -I${LOCAL_NAG_BLOCKS}
FC_DEBUG_FLAGS += -I${LOCAL_MKL_DIR} -I ${LOCAL_NAG_BLOCKS}

# Modules
MOD_OBJS = \
    ConstantsModule.o \
	ToolsModule.o \
	PropertyList.o \
	Ranking.o \
	NewTools.o \
	GaussianQuadrature.o \
	LinearAlgebra.o \
	GlobalModule.o \
	OutputModule.o \
	DataModule.o \
	LikelihoodModule.o

# Main program
PRG_OBJS = main.o

# Test modules
TEST_MOD_OBJS = tests/GaussianQuadratureTests.o tests/LinearAlgebraTests.o

# Test program
TEST_PRG_OBJS = tests/RunTests.o

OBJS = $(MOD_OBJS) $(PRG_OBJS)
DEBUG_OBJS = $(patsubst %.o, %_d.o, $(OBJS))
TEST_OBJS = $(TEST_MOD_OBJS) $(TEST_PRG_OBJS)

# Static pattern rules

.DEFAULT_GOAL := exe
$(OBJS) : %.o : %.f90
	$(FC) $(FC_FLAGS) -c $< -o $@

$(DEBUG_OBJS) : %_d.o : %.f90
	$(FC) $(FC_DEBUG_FLAGS) -c $< -o $@

$(TEST_OBJS) : %.o : %.f90
	$(FC) $(FC_FLAGS) -c $< -o $@

# Rules to build application

$(TGT).exe : $(OBJS)
	$(FC) $(FC_FLAGS) $^ -o $@ $(FL_LIBS)

$(TGT).dbg : $(DEBUG_OBJS)
	$(FC) $(FC_DEBUG_FLAGS) $^ -o $@ $(FL_LIBS)

RunTests.exe : $(MOD_OBJS) $(TEST_OBJS)
	$(FC) $(FC_FLAGS) $^ -o $@ $(FL_LIBS)

# Phony targets

.PHONY: exe dbg clean all

exe: $(TGT).exe

dbg: $(TGT).dbg

tests: RunTests.exe

clean:
	rm -f *.o *.mod *.dbg *.exe tests/*.o

all: clean dbg exe tests
