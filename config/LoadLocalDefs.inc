# LoadLocalDefs.inc - loads correct platform-specific definitions
#
#    SPARSE DEMAND
#
# modification history
# --------------------
# 2021MAR06 LN update paths
# 2020OCT18 LN update for myriad
#

# Select correct Platform-specific parameters
HOST_ID := $(shell hostname)
HOST_ID := $(if $(findstring login,$(HOST_ID)),myriad,$(HOST_ID))
USER    := $(shell whoami)

ifeq ($(HOST_ID),myriad)
  include $(WORK_DIR)/config/myriad.inc
else ifeq ($(HOST_ID),peacock.cs.ucl.ac.uk)
  include $(WORK_DIR)/config/cs.inc
else ifeq ($(HOST_ID),jake.local)
  include $(WORK_DIR)/config/jake.inc
else 
  include $(WORK_DIR)/config/local.inc
endif
