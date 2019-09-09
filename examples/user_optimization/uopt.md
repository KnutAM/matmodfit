# user optimization example for matmodfit

## How to build this example without CMake

### Windows: 
Create empty Visual Studio dll project named uopt_example, add source file uopt.f90. 
Ensure that 64-bit architecture is chosen and build
Copy the resulting dll file "uopt_example.dll" to a place in your path 
(Alternatively change input file to absolute path to this file without extension)

### Linux comp=ifort or gfortran
comp -shared -fPIC uopt.f90 -o uopt_example.so
Copy uopt_example.so to somewhere in LD_LIBRARY_PATH, or modify input file to give absolute path (w/o .so extension)

## Run example
From the folder containing uopt.inp and uopt_expdata.txt open cmd/terminal and type
matmodfit uopt.inp

