#!/bin/bash
# Script to build umats from repository MaterialModels

if ["$1$" == "" ]; then
    FC=ifort
else
    FC="$1"
fi;

# Go to directory outside the matmodfit git environment
cd ../..

# Download the MaterialModels repository
git clone git@github.com:KnutAM/MaterialModels.git

cd MaterialModels/scripts

python compile_general.py


# Add compiled location to bash_profile etc. 
cd ..
compiled_dir="$PWD/compiled"
cd scripts
./add_to_system_variable_if_not_already_there.sh LD_LIBRARY_PATH $compiled_dir