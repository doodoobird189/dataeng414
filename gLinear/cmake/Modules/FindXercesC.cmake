#
# Try to find XercesC by using the cmake find_package command, e.g.:
#
#   find_package(XercesC REQUIRED)
#
# For more details on find_package options, see:
#
#   cmake --help-command find_package
# 
# Once done, this will define
#
# XERCESC_FOUND - system has XercesC
# XERCESC_INCLUDE_DIRS - the XercesC include directories
# XERCESC_LIBRARIES - link these to use XercesC
#
# If XercesC cannot be found in the locations specified in this
# script, you can manually forward known locations for the XercesC
# library and header files, by defining the variables
# XERCESC_LIBRARY_SEARCH_PATH and XERCESC_INCLUDE_SEARCH_PATH
# respectively, prior to calling find_package, e.g.
#
# set(XERCESC_INCLUDE_SEARCH_PATH $ENV{HOME}/Developer/usr/include)
# set(XERCESC_LIBRARY_SEARCH_PATH $ENV{HOME}/Developer/usr/lib)
# find_package(XercesC REQUIRED)
#

if(XERCESC_INCLUDE_DIRS AND XERCESC_LIBRARIES)
  set(XERCESC_FIND_QUIETLY TRUE)
endif()

# include dir
find_path(
  XERCESC_INCLUDE_DIR
    sax2/ContentHandler.hpp
  PATHS
    # look in common install locations first
    /usr/local/include/xercesc
    /usr/include/xercesc
    /opt/include/xercesc
    ${CYGWIN_INSTALL_PATH}/include/xercesc
    # then in a possible user specified location
    ${XERCESC_INCLUDE_SEARCH_PATH}
  )

# finally the library itself
find_library(
  XERCESC_LIBRARY NAMES
    xerces-c
  PATHS
    # look in common install locations first
    /usr/local/lib
    /usr/lib
    /opt/lib
    ${CYGWIN_INSTALL_PATH}/lib
    # then in a possible user specified location
    ${XERCESC_LIBRARY_SEARCH_PATH}
  )

# handle the QUIETLY and REQUIRED arguments and set XERCESC_FOUND to
# TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  XercesC DEFAULT_MSG
  XERCESC_LIBRARY XERCESC_INCLUDE_DIR
  )

set(XERCESC_LIBRARIES ${XERCESC_LIBRARY})
set(XERCESC_INCLUDE_DIRS ${XERCESC_INCLUDE_DIR})

mark_as_advanced(XERCESC_LIBRARY XERCESC_INCLUDE_DIR)
