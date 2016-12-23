# LoadLocalDefs.inc - loads correct platform-specific definitions
#
# modification history
# --------------------
# 11JUN2016 LN  update to new paths on UCL and CS clusters
# 08DEC2014 LN  hostname doesn't work on cs cluster
# 15NOV2014 LPN edit to work on CS cluster
# 23JAN2014 LN  edit to work on legion
# 16jun2011 LPN changed hostname for IFS
# 08feb2008 bss fixed host specification logic
# 07feb2008 bss written.
#
# Define possible host ID's
UCL_HPC = 100
IFS_HPC = 200
CS_HPC = 300
LAPTOP = 400

# Select correct Platform-specific parameters
HOST_ID := $(shell hostname)
#HOST_ID = vic.local

ifeq ($(HOST_ID),drebin)
HOST_ID=$(IFS_HPC)
include ../config/config.ifs.hpc.inc
else ifeq ($(HOST_ID),bergerac)
HOST_ID=$(IFS_HPC)
include ../config/config.ifs.hpc.inc
else ifeq ($(HOST_ID),sherlock)
HOST_ID=$(IFS_HPC)
include ../config/config.ifs.hpc.inc
else ifeq ($(HOST_ID),wallander)
HOST_ID=$(IFS_HPC)
include ../config/config.ifs.hpc.inc
else ifeq ($(HOST_ID),wise.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),vic.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),jake.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),zeppo-22-1.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),zeppo-22-2.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),zeppo-22-3.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),zeppo-22-8.local)
HOST_ID=$(CS_HPC)
include /home/uctpln0/GeneralCode/config/config.cs.inc
else ifeq ($(HOST_ID),Lars-Nesheims-MacBook-Pro.local)
HOST_ID=$(LAPTOP)
else
HOST_ID=$(UCL_HPC)
include /data/uctpln0/GeneralCode/config/config.ucl.hpc.inc
endif
