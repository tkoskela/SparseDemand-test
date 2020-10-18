# LoadLocalDefs.inc - loads correct platform-specific definitions
#
#    SPARSE DEMAND
#
# modification history
# --------------------
# 2020OCT18 LN update for myriad
# 11JUN2016 LN  update to new paths on UCL and CS clusters
# 08DEC2014 LN  hostname doesn't work on cs cluster
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

GeneralCode=/home/uctpln0/GeneralCode
 
ifeq ($(HOST_ID),myriad)
  HOST_ID=$(MYRIAD)
  include $(GeneralCode)/config/myriad.inc 
else ifeq ($(HOST_ID,peacock.local)
  HOST_ID=$(CS_HPC)
  include $(GeneralCode)/config/cs.inc
else ifeq ($(HOST_ID),jake.local)
  HOST_ID=$(CS_HPC)
  include $(GeneralCode)/config/jake.inc
endif
