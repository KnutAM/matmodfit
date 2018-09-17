#!/bin/bash
# Script to build matmodfit and its material models
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

# Try to set the MKLROOT environment variable if not set,
# by running intel's mklvars.sh script
if [ -z "$MKLROOT" ]; then
	mklvars=$(find /opt/intel/* -name mklvars.sh)
	source $mklvars intel64
	export MKLROOT="$MKLROOT"
fi;


# Modify the ~./profile file to contain the necessary environmental variables
START_STRING="#Automatically added matmodfit directories"
LD_STRING='export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":'"$installdir"
PATH_STRING='export PATH="$PATH":'"$installdir"
END_STRING="#End of matmodfit directories"

# We need to check if .bash_profile exists. If it does, it will override content of .profile
# We should therefore not create it if it doesn't exist.
path_file="noadd"
if [ -f ~/.bash_profile ]; then
	echo "Do you want to add path entries to the ~/.bash_profile file? (yes/no)"
	read answer
	if [ "$answer" == "yes" ]; then
		echo "$START_STRING" >>~/.bash_profile
		echo "$LD_STRING" >>~/.bash_profile
		echo "$PATH_STRING" >>~/.bash_profile
		echo "$END_STRING" >>~/.bash_profile
		path_file="~/.bash_profile"
	fi;
else
	echo "Do you want to add path entries to the ~/.profile file? (yes/no)"
	read answer
	if [ "$answer" == "yes" ]; then
		echo "$START_STRING" >>~/.profile
		echo "$LD_STRING" >>~/.profile
		echo "$PATH_STRING" >>~/.profile
		echo "$END_STRING" >>~/.profile
		path_file="~/.profile"
	fi;
fi;



# Create installation directory
mkdir $installdir  # Create folder to contain all the finished software

# Build main program
rm -r build
mkdir build  # Create a build folder
cd build

FC=$compiler cmake ../src
make

cp matmodfit $installdir

echo "matmodfit has been built, building material models"
# Build material models
cd ../scripts   # Exit out of build folder
python bld_umats.py $installdir $compiler

echo "material models have been built, building examples"
python bld_examples.py $installdir $compiler


PATH="$installdir:$PATH"
LD_LIBRARY_PATH="$installdir:$LD_LIBRARY_PATH"

echo $PATH
echo $LD_LIBRARY_PATH

echo "Verification of installation"
matmodfit --version


#Script completed
echo "Setup completed"
if [ "$path_file" == "noadd" ]; then
	echo "Path not automatically added, you have to ensure this manually"
else
	echo "Type: source $path_file to add program to path directly"
	echo "This will happen automatically after next logout"
fi;
