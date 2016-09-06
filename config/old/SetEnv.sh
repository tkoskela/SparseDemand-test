# Revision history
# 2015SEP11 LN  update paths to intel 2015.1.133 compiler and to NAG Mark 25
#               also update NAG license server address
# add intel, openmpi and NAG to LD_LIBRARY_PATH
for x in /share/apps/intel/composer_xe_2015.1.133/compiler/lib/intel64 /share/apps/openmpi-1.10.0-intel-15.0.1/lib /share/apps/NAG/lib/fll6i25dcl/lib; do
  case ":$LD_LIBRARY_PATH:" in
    *":$x:"*) :;; #already there
    *) LD_LIBRARY_PATH="$x:$LD_LIBRARY_PATH";;
  esac
done

# add intel, openmpi and java to path
for x in /share/apps/intel/bin /share/apps/openmpi-1.10.0-intel-15.0.1/bin /share/apps/jdk1.8.0_25/bin; do
  case ":$PATH:" in
    *":$x:"*) :;; #already there
    *) PATH="$x:$PATH";;
  esac
done

# add intel man path to manpath
for x in /share/apps/intel/man/en_US; do
  case ":$MANPATH:" in
    *":$x:"*) :;; # already there
    *) MANPATH="$MANPATH:$x";;
  esac
done

export INTEL_LICENSE_FILE=28518@lic-intel.ucl.ac.uk
export NAG_KUSARI_FILE=10.1.1.1:
#export NAG_KUSARI_FILE=naglm-server.ucl.ac.uk:
export MANPATH
export LD_LIBRARY_PATH
export PATH

/share/apps/intel/bin/compilervars.sh intel64
