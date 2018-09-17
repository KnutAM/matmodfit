# matmodfit 
material model fitting program

## Download with git
"git clone git@github.com:KnutAM/matmodfit.git"

## Language version and compiler support
The Fortran 2008 standard is used
Compilers ifort (v. 16 or newer) or gfortran (v6 or newer)
Linux (gfortran or ifort)
Windows (ifort with Visual studio)

## Setup environment
Windows
* "mkdir build"
* "cd build"
* "cmake -A x64 ../src"

Linux
* "cd scripts"
* "chmod 744 install_matmodfit.sh"
* ./install_matmodfit.sh

## Coding standard
- Coding standard should be driven towards: http://www.fortran.com/Fortran_Style.pdf
- Lower case letters should primarily be used (exceptions are parameters)
- All modules should be private by default, and only export the necessary components. 
- All functions should have a short text directly under the function definition line, explaining the main purpose of the function
- Whenever a variable name is defined it should have a short explanation of what it is. 
- All double precision variables that are assigned a number, should be have d0 or equivalent appended to ensure that true double precision is invoked. 

### naming scheme
- Procedures: calculate_square(var_one, var_two) - Active names, describing their operation
- Variables: Normally like "this_is_a_variable_name_example"
- Parameters: Like variables, but with capital letters
- Modules: End with "_mod" ()
    - The file name should match module name 
    - If different files are used for different system configurations, append configuration identifier to the file name (e.g. _win)
- Types: End with "_typ"
- Pointers: End with "_ptr"

### Compiler warnings
All warnings and checks should be enabled during debug. 
Except for the cases described below, warnings should be treated as errors:
- "An array temporary was created" (ifort: -check:arg_temp_creation)
    If possible this warning should be avoided. However, it can be accepted (although still not preferred) if the following conditions are met:
    1) The variable/array is relatively small
    2) The number of calls are limited
    3) It is very difficult/cumbersome to avoid
    Note that this warning can occur if the wrong order of array indices are used, should be typical (:,k) and not (k,:) if passing a 1d array from a 2d. 
