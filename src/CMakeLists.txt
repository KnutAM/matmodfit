cmake_minimum_required (VERSION 3.0)
project (matmodfit)
enable_language (Fortran)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeScripts" CACHE STRING "Choose folder for cmake scripts")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# SYSTEM SPECIFICS
# Compiler flags needs to be added before introducing the add_executable/add_library
# Get compiler name
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
# Define logic variable for linux
if (UNIX AND NOT APPLE)
  set(LINUX TRUE)
else()
  set(LINUX FALSE)
endif()
# Define logic variable for architecture
if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(A64 TRUE)
else()
    set(A64 FALSE)
endif()

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("gfortran used")
    set(GFORT TRUE)
    set(IFORT FALSE)
    if (A64)
        add_compile_options(-m64 -ffree-line-length-none)
    else ()
        add_compile_options(-ffree-line-length-none)
    endif()
    add_compile_options("$<$<CONFIG:DEBUG>:-Wall>")
    add_compile_options("$<$<CONFIG:DEBUG>:-Wextra>")
    add_compile_options("$<$<CONFIG:DEBUG>:-Wno-tabs>")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    message ("ifort used")
    set(GFORT FALSE)
    set(IFORT TRUE)
    if (WIN32)
        add_compile_options("$<$<CONFIG:DEBUG>:/warn>")
        add_compile_options("$<$<CONFIG:DEBUG>:/check>")
    endif()
endif()

if (MSVC)
    add_compile_options(/names:lowercase /assume:underscore /fpp)
endif()


# == General library (gen) == #
# Note that some parts of this depends on compiler and os, see target_sources further down
add_library(gen STATIC 
	"gen/gen_util_mod.f90"
    "gen/types_mod.f90" 
    "gen/output_mod.f90"
    "gen/convert_mpar_mod.f90"
    "gen/nlopt_constants_mod.f90"              #nlopt constants
    "gen/umat_mod.f90"
	"gen/interface_mod.f90"
	"gen/usr_interface_mod.f90"
	"gen/load_dll_mod.f90"
    )
add_library(matmodfit::gen ALIAS gen)
#target_compile_options(gen PRIVATE ${CMAKE_DL_LIBS})

# == Input library (inp) == #
add_library(inp STATIC
    "inp/input_mod.f90"
    "inp/reading_utilities_mod.f90"
    "inp/sim_inp_mod.f90"
    "inp/opt_inp_mod.f90"
    "inp/check_input_mod.f90"
    )
add_library(matmodfit::inp ALIAS inp)
target_link_libraries(inp PRIVATE matmodfit::gen)

# == Simulation library (sim) == #
add_library(sim STATIC "${PROJECT_SOURCE_DIR}/sim/simulation_mod.f90")
add_library(matmodfit::sim ALIAS sim)
target_link_libraries(sim PRIVATE matmodfit::gen)
target_sources(sim
    PRIVATE
     "sim/sim_setup_mod.f90"
     "sim/sim_writing_util_mod.f90"
     "sim/sim_getincr_mod.f90"
     "sim/sim_error_mod.f90"
     "sim/sim_util_mod.f90"
     "sim/sim_tensor_mod.f90"
     "sim/atp/atp_mod.f90"
     "sim/atp/atp_util_mod.f90"
	 "sim/atp/atp_import_mod.f90"
     "sim/atp/atp_element_mod.f90"
     "sim/atp/atp_material_removal_mod.f90"
     "sim/mps/mps_mod.f90"
     "sim/mps/mps_util_mod.f90"
     )

# == Optimization library (opt) == #
add_library(opt STATIC 
     "opt/nlopt_objective.f90"
     "opt/nlopt_mod.f90"
	 "opt/optimization_mod.f90"
	 )
add_library(matmodfit::opt ALIAS opt)
target_link_libraries(opt PRIVATE matmodfit::gen)
target_link_libraries(opt PRIVATE matmodfit::sim)

# == Main executable (matmodfit) == #
add_executable(matmodfit "matmodfit.f90" "run_types_mod.f90")
target_link_libraries(matmodfit PRIVATE matmodfit::gen) # Not actually needed because it is given via all the others
target_link_libraries(matmodfit PRIVATE matmodfit::inp)
target_link_libraries(matmodfit PRIVATE matmodfit::opt)
target_link_libraries(matmodfit PRIVATE matmodfit::sim) # Not actually needed because it is given via matmodfit::opt
    

if (LINUX AND A64)
  message("Linux 64 bit detected")
  # Add special linux module
  target_sources(gen PRIVATE "${PROJECT_SOURCE_DIR}/gen/dll_mod_lin.f90")
  target_link_libraries(gen PRIVATE ${CMAKE_DL_LIBS})
  
  # Link to nlopt library
  # find_package(nlopt)
  find_library(nlopt_lib nlopt)
  target_link_libraries(opt PRIVATE "${nlopt_lib}")
#  set (NLOPTROOT "${PROJECT_SOURCE_DIR}/../lib" CACHE STRING "Please specify your nlopt path")
#  target_link_libraries(opt PRIVATE "${NLOPTROOT}/libnlopt.so")
  
  # Link to MKL library
  # Check if MKLROOT environment variable is set, otherwise try to set it using mklvars.sh
  set (MKLROOT $ENV{MKLROOT} CACHE STRING "Please specify your MKLROOT path")
  	
  target_include_directories(sim PRIVATE "${MKLROOT}/include")
  if (GFORT)
    target_link_libraries(sim PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
    target_sources(gen PRIVATE "${PROJECT_SOURCE_DIR}/gen/compiler_mod_gfort.f90")
  elseif (IFORT)
    target_link_libraries(sim PRIVATE "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
    target_sources(gen PRIVATE "${PROJECT_SOURCE_DIR}/gen/compiler_mod_ifort.f90")
  endif ()
  
elseif (WIN32) #Also true for 64 bit systems
  message("Windows detected")
  # Add special windows module
  target_sources(gen PRIVATE "${PROJECT_SOURCE_DIR}/gen/dll_mod_win.f90")
  
  if (IFORT)
    # Add special ifort module
    target_sources(gen PRIVATE "${PROJECT_SOURCE_DIR}/gen/compiler_mod_ifort.f90")
    if (NOT A64)
     # Link to mkl according to Intel's instructions
     target_link_libraries(sim PRIVATE "mkl_intel_c.lib mkl_sequential.lib mkl_core.lib")
     # Link to nlopt using the supplied binaries.
     set (NLOPTROOT "${PROJECT_SOURCE_DIR}/../nlopt_dll/x86" CACHE STRING "Please specify your nlopt path")
    elseif (A64)
     # Link to mkl according to Intel's instructions
     target_link_libraries(sim PRIVATE "mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib")
     # Link to nlopt using the supplied binaries.
     set (NLOPTROOT "${PROJECT_SOURCE_DIR}/../nlopt_dll/x64" CACHE STRING "Please specify your nlopt path")
    endif ()
  else()
    message("Only ifort compiler supported on Windows")
  endif()
  
  target_include_directories(opt PRIVATE ${NLOPTROOT})
  target_link_libraries(opt PRIVATE "${NLOPTROOT}/libnlopt-0.lib")
  
else()
  message("Only Windows and 64 bit Linux supported")
endif()
