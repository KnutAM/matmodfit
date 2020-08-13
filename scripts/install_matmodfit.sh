#!/bin/bash
# Script to build matmodfit
# Run from inside the script folder
# You may need to add permissions as chmod 744 install_matmodfit.sh first


# Setup default directories
cd ..	## Exit out of the script folder
installdir="$PWD/matmodfit"

# Ask user about which compiler to use
echo "Please specify which compiler to use: 1 ifort, 2 gfortran"
read answer
if [ "$answer" == "2" ]; then
    compiler="gfortran"
else
    compiler="ifort"
fi;

# Ask user for alternative installation directory
echo "Do you want to install in the default directory: $installdir ? (yes/no)"
read answer
if [ "$answer" == "no" ]; then
	echo "Please specify where you want to install:"
	read user_install_dir
	installdir=$user_install_dir
fi;

# Ask the user to install MaterialModels
echo "matmodfit require material routines to run, do you want to install the default from the MaterialModels repository (required for tests to work)? (yes/no)"
read install_MaterialModels


# Try to set the MKLROOT environment variable if not set,
# by running intel's mklvars.sh script
if [ -z "$MKLROOT" ]; then
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
	export MKLROOT="$MKLROOT"
fi;



if [[ :$PATH: == *:"$installdir":* ]] ; then
    PATH=$PATH:$installdir
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$installdir
fi;

chmod 744 add_to_system_variable_if_not_already_there.sh
./add_to_system_variable_if_not_already_there.sh PATH $installdir
./add_to_system_variable_if_not_already_there.sh LD_LIBRARY_PATH $installdir

# Create installation directory
mkdir $installdir  # Create folder to contain all the finished software

# Build main program
rm -r build
mkdir build  # Create a build folder
cd build

if [ -z "$nlopt_path" ]; then
    FC=$compiler cmake ../src
else
    FC=$compiler cmake -D nlopt_lib=$nlopt_path ../src
fi;

if [ $? -ne 0 ]; then
    echo "Could not setup build environment, check that nlopt and mkl was found"
    echo "If in a cluster environment you must load the required modules for mkl"
    echo "If nlopt is not found please specify the full path of libnlopt.so"
    echo "or press enter to quit:"
    read nlopt_path
    if [ $nlopt_path ]; then
        rm -rf *
        FC=$compiler cmake -D nlopt_lib=$nlopt_path ../src
        if [ $? -ne 0 ]; then
            echo "Setting up environment still failed, I don't know why it isn't working..."
            rm -r $installdir
            exit 1
        fi;
    else
        rm -r $installdir
        exit 1
    fi;
fi;

make
if [ $? -ne 0 ]; then
    echo "Building of matmodfit failed, please see error messages"
    rm -r $installdir
    exit 1
fi;

cp matmodfit $installdir

echo "matmodfit has been built, building examples"
# Build material models
cd ../scripts   # Exit out of build folder
python bld_examples.py $installdir $compiler

echo "examples built"
if [ "$install_MaterialModels" == "yes" ]; then
    echo "Installing material models"
    chmod 744 install_MaterialModels.sh
    ./install_MaterialModels.sh
fi;

echo "Verification of installation"
matmodfit --version

#Script completed
echo "Setup completed"