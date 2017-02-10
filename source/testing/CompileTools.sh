#!/bin/sh
#
#  script to compile TestTools.f90
#
# revision history
# 2015AUG14  LN  created.

#INTEL_MKL_DIR="/cm/shared/apps/intel/composer_xe_2015.1.133/mkl/include"
#NR_DIR="/data/uctpln0/GeneralCode/nr"
#TOOLS_DIR="/data/uctpln0/GeneralCode/tools"
#NAG_DIR="/cm/shared/apps/NAG/fll6i25dc"
INTEL_MKL_DIR="/share/apps/intel/composer_xe_2015.1.133/mkl/include"
NR_DIR="/home/uctpln0/GeneralCode/nr"
TOOLS_DIR="/home/uctpln0/GeneralCode/tools"
NAG_DIR="/share/apps/NAG/lib/fll6i25dcl"
NAG_MKL_DIR="${NAG_DIR}/mkl_intel64_11.2.0/lib"
LIBS="${NAG_DIR}/lib/libnag_mkl.a -Wl,--start-group ${NAG_MKL_DIR}/libmkl_core.a ${NAG_MKL_DIR}/libmkl_intel_lp64.a ${NAG_MKL_DIR}/libmkl_intel_thread.a -Wl,--end-group -liomp5 -lpthread -lm -shared-intel"
FC_FLAGS="-fpp -cxxlib -g -debug extended -parallel"
FC="ifort"
INCLUDES="-I${INTEL_MKL_DIR} -I${NAG_DIR}/nag_interface_blocks"
MOD1="${NR_DIR}/nrtype.f90 ${NR_DIR}/nrutil.f90 ${NR_DIR}/nr.f90 ${NR_DIR}/ran_state.f90 ${NR_DIR}/numerical_recipes.f90 ${TOOLS_DIR}/ToolsModule.f90" 
MOD2="../NewTools.f90"
INPUTS="TestTools.f90"

${FC} ${FC_FLAGS} ${INCLUDES} ${MOD1} ${MOD2} ${INPUTS} -o TestTools.dbg ${LIBS}
