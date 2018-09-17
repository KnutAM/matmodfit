# Scripts for matmodfit

## install_matmodfit.sh
Building and installation script for Linux. 
Will ask the user for:
    * Which compiler (ifort or gfortran)
    * Installation directory (install_dir)
    * If install_dir should be added to PATH/LD_LIBRARY_PATH via .bash_profile or .profile
    * If nlopt library isn't found, it will ask the user for this library path and retry to build. 
It first builds the main program before building the included material models and example user libraries. 
It uses bld_examples.py and bld_umats.py, hence python is required. 

## install_matmodfit.vbs
Copying of files on Windows. Must be put in the same folder as the compiled matmodfit.exe program and the other dll's. 
Will put generate a start menu shortcut to the matmodfit.bat script that starts a cmd with the folder containing install_matmodfit.vbs in path. 
Will also copy the manual to the start menu.

## bld_examples.py and bld_umats.py
Builds the examples or included material models. 
First argument is output directory (default is ../build_[examples/umats])
Second argument is compiler (default is ifort)

## get_mklroot.sh
Will try to get the MKL environment variables and print the $MKLROOT

##matmodfit.bat
Start a cmd window with the path to the matmodfit program in path (for users without rights to add locations to path)