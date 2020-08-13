#!/bin/bash
cd ../..    # Go outside matmodfit folder

# Load modules
module load intel CMake Python

# Install nlopt
echo "Installing nlopt"
git clone git://github.com/stevengj/nlopt
cd nlopt
mkdir build
cd build
CC=icc CXX=icpc FC=ifort cmake -DCMAKE_INSTALL_PREFIX="$(pwd)/../install" ..
make
make install
cd ../..
export nlopt_path="$(pwd)/nlopt/install/lib64/libnlopt.so"

# Install matmodfit
echo "Installing matmodfit"
# git clone https://github.com/KnutAM/matmodfit.git (assume already downloaded)
cd matmodfit/scripts
chmod 744 install_matmodfit.sh

# Set the replies to avoid user interaction
echo "1" > replies.txt      # compiler=ifort (Replace for first to ensure empty file)
echo "yes" >> replies.txt   # default directory (Append for the remaining answers)
echo "yes" >> replies.txt   # install MaterialModels (Append for the remaining answers)

# Run the installation script
./install_matmodfit.sh < replies.txt    # Use the lines in replies.txt as answers. 

rm replies.txt # Remove the replies.txt file to avoid clutter


if [ -f "$HOME/.bash_profile" ]
then
    echo "Running test cases"
    source "$HOME/.bash_profile"
    cd ../testing
    python run_test_cases.py
    python clean_test.py
fi