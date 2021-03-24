# LoadLocalDefs.inc - loads correct platform-specific definitions
#
#    SPARSE DEMAND
#
# modification history
# --------------------
# 2021MAR06 LN update paths
# 2020OCT18 LN update for myriad
#
# Define possible host ID's
ECON_HPC = 100
IFS_HPC = 200
CS_HPC = 300
LAPTOP = 400
MYRIAD = 500

# Select correct Platform-specific parameters
HOST_ID := $(shell hostname)
HOST_ID := $(if $(findstring login,$(HOST_ID)),myriad,$(HOST_ID))
USER    := $(shell whoami)

ifeq ($(HOST_ID),myriad)
  HOST_ID=$(MYRIAD)
  include $(WorkDir)/config/myriad.inc
else ifeq ($(HOST_ID,peacock.cs.ucl.ac.uk)
  HOST_ID=$(CS_HPC)
  include $(WorkDir)/config/cs.inc
else ifeq ($(HOST_ID),jake.local)
  HOST_ID=$(CS_HPC)
  include $(WorkDir)/config/jake.inc
endif
