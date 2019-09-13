#!/bin/bash
# Note: This script is an auxillary script and is not used by install_matmodfit (this code is written directly in that file)
mklvars=$(find /opt/intel/* -name mklvars.sh)
if [ -z "$mklvars" ]
then
    mklvars=$(find $HOME/intel/* -name mklvars.sh)
fi
if [ -z "$mklvars" ]
then
    echo "Cannot detect that mkl is installed, please give the path to the ''intel/'' folder containing the mkl libraries"
    echo "E.g. /opt/intel/"
    read user_intel_dir
    mklvars=$(find $user_intel_dir -name mklvars.sh)
fi

if [ -z "$mklvars" ]
then
    echo "Still could not detect mkl in the specified folder, the installation will fail. Please verify that mklvars.sh is located in some subfolder of the specified intel folder"
    echo "Press enter to quit"
    read 
    return 1
fi

source $mklvars intel64
echo "MKL found"
echo $MKLROOT

