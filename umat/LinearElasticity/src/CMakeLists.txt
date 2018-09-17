cmake_minimum_required (VERSION 3.0)
project (LinearElasticExample)
enable_language (Fortran)

## This CMakeLists.txt file only creates an example shared library umat - Linear elasticity
## To build for standard system on your computer, call (from build folder): cmake ../src
## However, it is important that if you build (or download) mpfit with 64 bit, then the umats must also be in 64 bit
## To build for 64 bit and visual studio, call (from build folder): cmake -G "Visual Studio 11 2012 Win64" ../src
## if you have this particular version of Visual studio installed.

# Set some compile options needed to compile LinearElastic.for with gfortran
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	add_compile_options(-ffree-line-length-none -ffree-form)
endif()

# Add files
add_library(LinearElastic SHARED "LinearElastic.for")
