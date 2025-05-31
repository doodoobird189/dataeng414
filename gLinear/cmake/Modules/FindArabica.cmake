#
# Try to find Arabica by using the cmake find_package command, e.g.:
#
#   find_package(Arabica REQUIRED)
#
# For more details on find_package options, see:
#
#   cmake --help-command find_package
# 
# Once done, this will define
#
# ARABICA_FOUND - system has Arabica
# ARABICA_INCLUDE_DIRS - the Arabica include directories
# ARABICA_LIBRARIES - link these to use Arabica
#
# If Arabica cannot be found in the locations specified in this
# script, you can manually forward known locations for the Arabica
# library and header files, by defining the variables
# ARABICA_LIBRARY_SEARCH_PATH and ARABICA_INCLUDE_SEARCH_PATH
# respectively, prior to calling find_package, e.g.
#
# set(ARABICA_INCLUDE_SEARCH_PATH $ENV{HOME}/Developer/usr/include)
# set(ARABICA_LIBRARY_SEARCH_PATH $ENV{HOME}/Developer/usr/lib)
# find_package(Arabica REQUIRED)
#

if(ARABICA_INCLUDE_DIRS AND ARABICA_LIBRARIES)
  set(ARABICA_FIND_QUIETLY TRUE)
endif()

# include dir
find_path(
  ARABICA_INCLUDE_DIR
    SAX/ContentHandler.hpp
  PATHS
    # look in common install locations first
    /usr/local/include/arabica
    /usr/include/arabica
    /opt/include/arabica
    ${CYGWIN_INSTALL_PATH}/include/arabica
    # then in a possible user specified location
    ${ARABICA_INCLUDE_SEARCH_PATH}
  )

# finally the library itself
find_library(
  ARABICA_LIBRARY NAMES
    arabica
  PATHS
    # look in common install locations first
    /usr/local/lib
    /usr/lib
    /opt/lib
    ${CYGWIN_INSTALL_PATH}/lib
    # then in a possible user specified location
    ${ARABICA_LIBRARY_SEARCH_PATH}
  )

# handle the QUIETLY and REQUIRED arguments and set ARABICA_FOUND to
# TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Arabica DEFAULT_MSG
  ARABICA_LIBRARY ARABICA_INCLUDE_DIR
  )

set(ARABICA_LIBRARIES ${ARABICA_LIBRARY})
set(ARABICA_INCLUDE_DIRS ${ARABICA_INCLUDE_DIR})

mark_as_advanced(ARABICA_LIBRARY ARABICA_INCLUDE_DIR)
