#
# Try to find SPTK by using the cmake find_package command, e.g.:
#
#   find_package(SPTK REQUIRED)
#
# For more details on find_package options, see:
#
#   cmake --help-command find_package
# 
# Once done, this will define
#
# SPTK_FOUND        - system has SPTK
# SPTK_INCLUDE_DIRS - the SPTK include directories
# SPTK_LIBRARIES    - link these to use SPTK
#
# If SPTK cannot be found in the locations specified in this script,
# you can manually forward known locations for the SPTK libraries,
# header files and executables, by defining the variables
# SPTK_LIBRARY_SEARCH_PATH and SPTK_INCLUDE_SEARCH_PATH, respectively,
# prior to calling find_package, e.g.
#
# set(SPTK_INCLUDE_SEARCH_PATH $ENV{HOME}/Developer/SPTK/include)
# set(SPTK_LIBRARY_SEARCH_PATH $ENV{HOME}/Developer/SPTK/lib)
# find_package(SPTK REQUIRED)
#
# It is highly recommended that you install SPTK to the default location, i.e.
# /usr/local/SPTK
#
# installing SPTK
# cd /tmp
# wget http://sourceforge.net/projects/sp-tk/files/SPTK/SPTK-3.6/SPTK-3.6.tar.gz
# cd SPTK-3.6
# ./configure
# make
# make install OR sudo make install
#

if(SPTK_INCLUDE_DIRS AND SPTK_LIBRARIES)
  set(SPTK_FIND_QUIETLY TRUE)
endif()

find_path(
  SPTK_INCLUDE_DIR
    SPTK.h
  PATHS
    /usr/local/SPTK/include
    /opt/SPTK/include
    ${CYGWIN_INSTALL_PATH}/SPTK/include
    ${SPTK_INCLUDE_SEARCH_PATH}
    ../SPTK-3.6/include
    ../SPTK-3.5/include
  )

find_library(
  SPTK_LIBRARY NAMES
    SPTK
  PATHS
    /usr/local/SPTK/lib
    /opt/SPTK/lib
    ${CYGWIN_INSTALL_PATH}/SPTK/lib
    ${SPTK_INCLUDE_SEARCH_PATH}
    ../SPTK-3.6/lib
    ../SPTK-3.5/lib
  )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
  SPTK DEFAULT_MSG
  SPTK_LIBRARY SPTK_INCLUDE_DIR
  )

set(SPTK_LIBRARIES ${SPTK_LIBRARY})
set(SPTK_INCLUDE_DIRS ${SPTK_INCLUDE_DIR})

mark_as_advanced(SPTK_INCLUDE_DIR SPTK_LIBRARY)

# set paths in which to search for SPTK executables
set(SPTK_POSSIBLE_BIN_PATHS
  # look in common install locations first
  /usr/local/SPTK/bin
  /opt/SPTK/bin
  ${CYGWIN_INSTALL_PATH}/SPTK/bin
  ${SPTK_BIN_SEARCH_PATH}
  ../SPTK-3.6/bin
  ../SPTK-3.5/bin
  )

# check that we have all the required SPTK executables
if(SPTK_FIND_COMPONENTS)
  foreach(component ${SPTK_FIND_COMPONENTS})
    string(TOUPPER ${component} _COMPONENT)
    find_program(${_COMPONENT}
      NAMES ${component}
      PATHS ${SPTK_POSSIBLE_BIN_PATHS}
      NO_DEFAULT_PATH
      )
    if(${_COMPONENT})
      find_package_handle_standard_args(SPTK DEFAULT_MSG ${_COMPONENT})
    else()
      message(FATAL_ERROR "SPTK: Can't find ${component}")
    endif()
  endforeach()
endif()
