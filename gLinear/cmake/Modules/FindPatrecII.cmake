#
# Try to find patrecII by using the cmake find_package command, e.g.:
#
#   find_package(patrecII REQUIRED)
#
# For more details on find_package options, see:
#
#   cmake --help-command find_package
#
# Once done, this will define
#
# PATRECII_FOUND - system has patrecII
# PATRECII_INCLUDE_DIRS - the patrecII include directories
# PATRECII_LIBRARIES - link these to use patrecII
# PATRECII_SCRIPT_PATH - the path containing patrecII scripts
#
# If patrecII cannot be found in the locations specified in this
# script, you can manually forward known locations for the patrecII
# library and header files, by defining the variables
# PATRECII_LIBRARY_SEARCH_PATH and PATRECII_INCLUDE_SEARCH_PATH
# respectively, prior to calling find_package, e.g.
#
# set(PATRECII_INCLUDE_SEARCH_PATH $ENV{HOME}/Developer/usr/include)
# set(PATRECII_LIBRARY_SEARCH_PATH $ENV{HOME}/Developer/usr/lib)
# find_package(patrecII REQUIRED)
#

if(PATRECII_INCLUDE_DIRS AND PATRECII_LIBRARIES)
  set(PATRECII_FIND_QUIETLY TRUE)
endif()

# include dir
find_path(
  PATRECII_INCLUDE_DIR
    _PATRECII_SOURCE_ROOT_
  PATHS
    # look in common install locations first
    /usr/local/include
    /usr/include
    /opt/include
    ${CYGWIN_INSTALL_PATH}/include
    # then in a possible user specified location
    ${PATRECII_INCLUDE_SEARCH_PATH}
    # or a relative location
    ../patrecII/src
    # or in documented/standard location
    $ENV{HOME}/devel/patrecII/src
  )

# script dir
find_path(
  PATRECII_SCRIPT_DIR
    patids.sh
  PATHS
    # look in common dev locations
    ../patrecII/scripts
    $ENV{HOME}/devel/patrecII/scripts
  )

# finally the library itself
find_library(
  PATRECII_LIBRARY NAMES
    patrecII
  PATHS
    # look in common install locations first
    /usr/local/lib
    /usr/lib
    /opt/lib
    ${CYGWIN_INSTALL_PATH}/lib
    # then in a possible user specified location
    ${PATRECII_LIBRARY_SEARCH_PATH}
    # or relative locations
    ../patrecII/obj.i386/lib
    ../patrecII/build/${TARGET_PLATFORM}/lib
    # or in documented/standard locations
    $ENV{HOME}/devel/patrecII/obj.i386/lib
    $ENV{HOME}/devel/patrecII/build/${TARGET_PLATFORM}/lib
  )

# handle the QUIETLY and REQUIRED arguments and set PATRECII_FOUND to
# TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  patrecII DEFAULT_MSG
  PATRECII_LIBRARY PATRECII_INCLUDE_DIR
  )

set(PATRECII_LIBRARIES ${PATRECII_LIBRARY})
set(PATRECII_INCLUDE_DIRS ${PATRECII_INCLUDE_DIR})
set(PATRECII_SCRIPT_PATH ${PATRECII_SCRIPT_DIR})

mark_as_advanced(PATRECII_LIBRARY PATRECII_INCLUDE_DIR)

# set paths in which to search for patrecII executables
set(PatrecII_POSSIBLE_BIN_PATHS
  # look in common install locations first
  /usr/bin
  /usr/local/bin
  /opt/bin
  ${CYGWIN_INSTALL_PATH}/bin
  # then in a possible user specified location
  ${PATRECII_BINARY_SEARCH_PATH}
  # or relative locations
  ../patrecII/obj.i386/bin
  ../patrecII/build/${TARGET_PLATFORM}/bin
  # or in documented/standard locations
  $ENV{HOME}/devel/patrecII/obj.i386/bin
  $ENV{HOME}/devel/patrecII/build/${TARGET_PLATFORM}/bin
  )

include(FindPackageHandleStandardArgs)

# check that we have all the required patrecII executables
if(PatrecII_FIND_COMPONENTS)
  foreach(component ${PatrecII_FIND_COMPONENTS})
    string(TOUPPER ${component} _COMPONENT)
    find_program(${_COMPONENT}
      NAMES ${component}
      PATHS ${PatrecII_POSSIBLE_BIN_PATHS}
      )
    if(${_COMPONENT})
      find_package_handle_standard_args(PatrecII DEFAULT_MSG ${_COMPONENT})
    else()
      message(FATAL_ERROR "PatrecII: Can't find ${component}")
    endif()
  endforeach()
endif()

include(CMakeParseArguments)

# call this function to create a custom patrecII command
#
# COMMAND       the patrec executable (required)
# DEPENDS       file and target dependencies
# OUTPUT        file(s) to be produced as a result of running patrec command
# LIST_OUTPUT   specify rules to generate a list file after running patrec command
# ARGS          the argument sequence that will be piped to the patrec command
# DESCRIPTION   a string that will be displayed while running this command
#
# must specify either an OUTPUT or a LIST_OUTPUT
#
# NB: calling this function will not actually run the command. the command
# will only run when one of its outputs is made a dependency of a target.
# Refer to CMake's documentation for more detail on custom commands and targets
#
# example:
#
# patrec_add_custom_command(
#   COMMAND objfun
#   OUTPUT  segio.t02.gz
#   ARGS    2 TRC02 * - 1 segio.t02.gz 0
#   )
#
# add_custom_target(
#   segio ALL
#   DEPENDS segio.t02.gz
#   )
#
function(patrec_add_custom_command)
  set(options )
  set(oneValueArgs COMMAND DESCRIPTION)
  set(multiValueArgs DEPENDS OUTPUT LIST_OUTPUT ARGS)

  cmake_parse_arguments(PACC "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(PACC_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown keywords given to PATREC_ADD_CUSTOM_COMMAND(): \"${PACC_UNPARSED_ARGUMENTS}\"")
  endif()

  if(NOT PACC_FAIL_MESSAGE)
    set(PACC_FAIL_MESSAGE  "DEFAULT_MSG")
  endif()

  if(NOT PACC_COMMAND)
    message(FATAL_ERROR "Please specify COMMAND to patrec_add_custom_command!")
  endif()

  string(TOUPPER ${PACC_COMMAND} COMPONENT)
  set(PACC_COMMAND ${${COMPONENT}})
  if(NOT PACC_COMMAND)
    message(FATAL_ERROR "${COMPONENT} not found. Remember to add to list of components in FIND_PACKAGE statement")
  endif()

  string(REPLACE ";" "-" PACC_LOG "${PACC_OUTPUT}")
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/logs)
  set(PACC_LOG "${CMAKE_BINARY_DIR}/logs/${PACC_LOG}.log")

  string(REPLACE ";" " " PACC_ARGS "${PACC_ARGS}")

  if(PACC_LIST_OUTPUT)
    list(LENGTH PACC_LIST_OUTPUT LIST_OUTPUT_LEN)
    if(NOT LIST_OUTPUT_LEN EQUAL 3)
      message(FATAL_ERROR "LIST_OUTPUT received ${LIST_OUTPUT_LEN} parameters: ${PACC_LIST_OUTPUT}, but 3 parameters required.")
    endif()
    list(GET PACC_LIST_OUTPUT 0 LIST_FILE)
    list(GET PACC_LIST_OUTPUT 1 LIST_PATH)
    list(GET PACC_LIST_OUTPUT 2 LIST_PATTERN)
    get_filename_component(PACC_LOG ${LIST_FILE} NAME)
    set(PACC_LOG "${CMAKE_BINARY_DIR}/logs/${PACC_LOG}.log")
    if(NOT PACC_DESCRIPTION)
      set(PACC_DESCRIPTION "Generating ${LIST_FILE}")
    endif()
    set(OUTPUT_DIRECTION >> ${PACC_LOG})

    add_custom_command(
      OUTPUT ${PACC_OUTPUT} ${LIST_FILE}
      DEPENDS ${PACC_DEPENDS}
      COMMAND date "+======================== %F %T ========================" ${OUTPUT_DIRECTION}
      COMMAND echo "echo ${PACC_ARGS} | ${PACC_COMMAND}" ${OUTPUT_DIRECTION}
      COMMAND echo "=====================================================================" ${OUTPUT_DIRECTION}
      COMMAND echo ${PACC_ARGS} | ${PACC_COMMAND} ${OUTPUT_DIRECTION} 2>&1
      COMMAND ${FIND} ${LIST_PATH} -name ${LIST_PATTERN} | sort > ${LIST_FILE}
      COMMAND ${BADLST} ${LIST_FILE}
      COMMENT "${PACC_DESCRIPTION} ... \n$ echo ${PACC_ARGS} | ${PACC_COMMAND} ${OUTPUT_DIRECTION}"
      VERBATIM
      )
  else()
    if(NOT PACC_OUTPUT)
      message(FATAL_ERROR "No OUTPUT nor LIST_OUTPUT specified!")
    endif()
    list(GET PACC_OUTPUT 0 PACC_LOG)
    get_filename_component(PACC_LOG ${PACC_LOG} NAME)
    set(PACC_LOG "${CMAKE_BINARY_DIR}/logs/${PACC_LOG}.log")
    if(NOT PACC_DESCRIPTION)
      set(PACC_DESCRIPTION "Generating ${PACC_OUTPUT}")
    endif()
    set(OUTPUT_DIRECTION >> ${PACC_LOG})

    add_custom_command(
      OUTPUT ${PACC_OUTPUT}
      DEPENDS ${PACC_DEPENDS}
      COMMAND date "+======================== %F %T ========================" ${OUTPUT_DIRECTION}
      COMMAND echo "echo ${PACC_ARGS} | ${PACC_COMMAND}" ${OUTPUT_DIRECTION}
      COMMAND echo "=====================================================================" ${OUTPUT_DIRECTION}
      COMMAND echo ${PACC_ARGS} | ${PACC_COMMAND} ${OUTPUT_DIRECTION} 2>&1
      COMMENT "${PACC_DESCRIPTION} ... \n$ echo ${PACC_ARGS} | ${PACC_COMMAND} ${OUTPUT_DIRECTION}"
      VERBATIM
      )
  endif()
endfunction()
