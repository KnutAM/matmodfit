# - Try to find NLOPT
# Based on: https://gitlab.kitware.com/cmake/community/wikis/doc/tutorials/How-To-Find-Libraries
# See also: https://cmake.org/cmake/help/v3.12/manual/cmake-developer.7.html?#a-sample-find-module
# Once done this will define
#  NLOPT_FOUND - System has NLOPT
#  NLOPT_INCLUDE_DIRS - The NLOPT include directories
#  NLOPT_LIBRARIES - The libraries needed to use NLOPT
#  NLOPT_DEFINITIONS - Compiler switches required for using NLOPT

find_package(PkgConfig)
pkg_check_modules(PC_NLOPT QUIET nlopt)
set(NLOPT_DEFINITIONS ${PC_NLOPT_CFLAGS_OTHER})

find_path(NLOPT_INCLUDE_DIR nlopt.h
          HINTS ${PC_NLOPT_INCLUDEDIR} ${PC_NLOPT_INCLUDE_DIRS} $ENV{LD_LIBRARY_PATH} $HOME/install
          )

find_library(NLOPT_LIBRARY NAMES libnlopt nlopt
             HINTS ${PC_NLOPT_LIBDIR} ${PC_NLOPT_LIBRARY_DIRS} $HOME/install)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NLOPT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(nlopt  DEFAULT_MSG
                                  NLOPT_LIBRARY NLOPT_INCLUDE_DIR)

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY )

set(NLOPT_LIBRARIES ${NLOPT_LIBRARY} )
set(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR} )
