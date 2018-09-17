# user simulation example for matmodfit

## How to build this example
### Windows: 
Create empty Visual Studio dll project named usim_example, add source file usim.f90. 
Ensure that 64-bit architecture is chosen and build
Copy the resulting dll file "usim_example.dll" to a place in your path 
(Alternatively change input file to absolute path to this file without extension)

### Linux comp=ifort or gfortran
comp -shared -fPIC usim.f90 -o usim_example.so
Copy usim_example.so to somewhere in LD_LIBRARY_PATH, or modify input file to give absolute path (w/o .so extension)

## Run example
From the folder containing usim.inp in cmd/terminal type
matmodfit usim.inp
