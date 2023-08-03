# Bash shell commands to define IDL environment variables and aliases.
. /Applications/harris/idl/bin/idl_setup.bash

# start up IDL path
export IDL_PATH=+$IDL_DIR/examples:+$IDL_DIR/lib

# local user's IDL folder
export IDL=$HOME/idl

# NIRC2_IDL - NIRC2 data reduction & analysis
export NIRC2_DIR=$IDL/pipelines/nirc2
export IDL_PATH=${IDL_PATH}:+$NIRC2_DIR/pro

# Astro, Coyote, MPFIT packages
#export IDL_PATH=${IDL_PATH}:+$IDL/astron/pro:+$IDL/coyote:+$IDL/mpfit

# IDLUTILS (includes Astrolib, Coyote, & MPFIT)
export IDLUTILS_DIR=$IDL/idlutils
export IDL_PATH=${IDL_PATH}:+$IDLUTILS_DIR/goddard/pro:+$IDLUTILS_DIR/pro
export PATH=$IDLUTILS_DIR/bin:$PATH

# Other useful IDL programs
export IDL_PATH=${IDL_PATH}:+$NIRC2_DIR/utilities
