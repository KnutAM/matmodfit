cmake_minimum_required (VERSION 3.0)
project (StVernantExample)
enable_language (Fortran)

## This CMakeLists.txt file only creates an example shared library umat - StVernant type elasticity
## To build for 64 bit and visual studio, call (from build folder): cmake -G "Visual Studio 11 2012 Win64"
## if you have this particular version of Visual studio installed.
#
## MKL
## To show linking to mkl, this is included. However, if mkl is not installed on your system,
## or you have troubles finding it, this can be resolved by commenting out the dgemm call in StVernant.for
## in the function v9x9_dd_v9x9, and replace with the suggested matmul function call.

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
endif()


# Add files
add_library(Vernant SHARED "StVernant.for")

# Link to MKL library
if (MSVC)   #If Microsoft Visual compiler
    if (A64)
        target_link_libraries(Vernant "mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib")
    else()
        target_link_libraries(Vernant "mkl_intel_c.lib mkl_sequential.lib mkl_core.lib")
    endif()
elseif (A64)
    if (UNIX AND NOT APPLE) # I.e. linux
        set (MKLROOT "/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl" CACHE STRING "Please specify your MKLROOT path")
        target_include_directories(Vernant PRIVATE "${MKLROOT}/include")
        if (GFORT)
            target_link_libraries(Vernant PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        elseif (IFORT)
            target_link_libraries(Vernant PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
        endif ()
    endif()
else()
    message("System configuration not supported")
endif()
