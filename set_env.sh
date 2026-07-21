#!/bin/bash

# 1. Initialize module command (ensures script works in all shell modes)
source /etc/profile.d/modules.sh

echo ">>> Setting up Environment for CLAS12 RG-D EXTERNALS <<<"

# 2. Use JLab standard module path
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles

# 3. Load CLAS12 module
# This automatically loads cernlib/2023, ROOT, and other dependencies
module load clas12

# 4. Set Project Variables
export EXTERNDIR=$(pwd)

# 5. Fix for Makefile: Define CERN_ROOT
# The module sets $CERN and $CERN_LEVEL, but the Makefile expects $CERN_ROOT.
if [ -n "$CERN" ] && [ -n "$CERN_LEVEL" ]; then
    export CERN_ROOT=${CERN}/${CERN_LEVEL}
else
    echo "WARNING: CERN variables not detected. Make sure 'clas12' loaded correctly."
fi

echo ""
echo " > Environment Ready."
echo " > EXTERNDIR: $EXTERNDIR"
echo " > CERN_ROOT: $CERN_ROOT"