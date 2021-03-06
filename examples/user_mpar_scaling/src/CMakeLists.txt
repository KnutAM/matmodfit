# Setup with visual studio: <cmake -A x64 ../src>   (To build 64 bit)
# Setup with make (linux):  <cmake ../src>
# Build from command line: <cmake --build . --config RELEASE>
cmake_minimum_required (VERSION 3.0)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
        
project(umpscale)
enable_language (Fortran)

# Get system information
if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(A64 TRUE)   # 64 bit compilation system - needs to be same as for mpsim
else()
    set(A64 FALSE)  # 32 bit compilation system - needs to be same as for mpsim
endif()
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set(GFORT TRUE)
    set(IFORT FALSE)
    if (A64)
        add_compile_options(-m64 -ffree-line-length-none -ffree-form)
    else ()
        add_compile_options(-ffree-line-length-none -ffree-form)
    endif()
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    set(GFORT FALSE)
    set(IFORT TRUE)
	add_compile_options("$<$<CONFIG:DEBUG>:/warn>")
    add_compile_options("$<$<CONFIG:DEBUG>:/check>")
endif()


# Add files
add_library(umpscale SHARED "chab_scale.f90")
    
set_target_properties(umpscale PROPERTIES OUTPUT_NAME "user_mpar_scale_example")
set_target_properties(umpscale PROPERTIES PREFIX "")
