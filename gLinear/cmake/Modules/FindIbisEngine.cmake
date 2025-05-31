#
# Try to find IbisEngine by using the cmake find_package command, e.g.:
#
#   find_package(IbisEngine COMPONENTS
#     flite_fiona2 flite_cw_eng_za flite HTSEngine
#   )
#
# For more details on find_package options, see:
#
#   cmake --help-command find_package
# 
# Once done, this will define
#
# IBISENGINE_FOUND        - system has IbisEngine
# IBISENGINE_INCLUDE_DIRS - the Ibis include directories
# IBISENGINE_LIBRARIES    - link these to use IbisEngine. includes flite, HTSEngine,
#                           language packs and voices
# IBISENGINE_EXECUTABLES  - IbisEngine executables flite and hts_engine defined also as
#                           FLITE and ENGINE respecively (HTS voice building use ENGINE)
#
# If IbisEngine cannot be found in the locations specified in this
# script, you can manually forward known locations for the IbisEngine
# libraries, header files and executables, by defining the variables
# IBISENGINE_LIBRARY_SEARCH_PATH, IBISENGINE_INCLUDE_SEARCH_PATH and
# IBISENGINE_EXECUTABLE_SEARCH_PATH respectively, prior to calling
# find_package, e.g.
#
# set(IBISENGINE_INCLUDE_SEARCH_PATH $ENV{HOME}/Developer/usr/include)
# set(IBISENGINE_LIBRARY_SEARCH_PATH $ENV{HOME}/Developer/usr/lib)
# set(IBISENGINE_EXECUTABLE_SEARCH_PATH $ENV{HOME}/Developer/usr/bin)
# find_package(IbisEngine COMPONENTS flite_fiona2 flite_cw_eng_za flite HTSEngine)
#

if(IBISENGINE_INCLUDE_DIRS AND IBISENGINE_LIBRARIES)
  set(IBISENGINE_FIND_QUIETLY TRUE)
endif()

# find include directories
foreach(header flite.h cw_flite.h HTS_engine.h)
  find_path(
    IBISENGINE_${header}_DIR
    NAMES
      ${header}
    PATHS
      # look in common install locations first
      /usr/local/include
      /usr/include
      /opt/include
      ${CYGWIN_INSTALL_PATH}/include
      # then in a possible user specified location
      ${IBISENGINE_INCLUDE_SEARCH_PATH}
      # or relative locations
      ../ibis_engine/src/cw_flite/include
      ../ibis_engine/src/3rdparty/flite/include
      ../ibis_engine/src/3rdparty/hts_engine/include
      # or in documented/standard locations
      $ENV{HOME}/devel/ibis_engine/src/cw_flite/include
      $ENV{HOME}/devel/ibis_engine/src/3rdparty/flite/include
      $ENV{HOME}/devel/ibis_engine/src/3rdparty/hts_engine/include
    )
  mark_as_advanced(IBISENGINE_${header}_DIR)
  list(APPEND IBISENGINE_INCLUDE_DIR ${IBISENGINE_${header}_DIR})
endforeach()

# find libraries
# find_package(Ibis COMPONENTS flite_fiona2 flite_cw_eng_za flite HTSEngine REQUIRED)
if(IbisEngine_FIND_COMPONENTS)
  foreach(component ${IbisEngine_FIND_COMPONENTS})
    find_library(
      IBISENGINE_${component}_LIBRARY
      NAMES
        ${component}
      PATHS
        # look in common install locations first
        /usr/local/lib
        /usr/lib
        /opt/lib
        ${CYGWIN_INSTALL_PATH}/lib
        # then in a possible user specified location
        ${IBISENGINE_LIBRARY_SEARCH_PATH}
        # or relative locations
        ../ibis_engine/build/${TARGET_PLATFORM}/lib
        # or in documented/standard locations
        $ENV{HOME}/devel/ibis_engine/build/${TARGET_PLATFORM}/lib
      )
    mark_as_advanced(IBISENGINE_${component}_LIBRARY)
    list(APPEND IBISENGINE_LIBRARY ${IBISENGINE_${component}_LIBRARY})
  endforeach()
endif()

# find executables
foreach(executable hts_engine flite)
  find_program(
    IBISENGINE_${executable}_EXECUTABLE
    NAMES
      ${executable}
    PATHS
      # look in common install locations first
      /usr/local/bin
      /usr/bin
      /opt/bin
      ${CYGWIN_INSTALL_PATH}/bin
      # then in a possible user specified location
      ${IBISENGINE_EXECUTABLE_SEARCH_PATH}
      # or relative location
      ../ibis_engine/build/${TARGET_PLATFORM}/bin
      # or in documented/standard location
      $ENV{HOME}/devel/ibis_engine/build/${TARGET_PLATFORM}/bin
    )
  mark_as_advanced(IBISENGINE_${executable}_EXECUTABLE)
  list(APPEND IBISENGINE_EXECUTABLE ${IBISENGINE_${executable}_EXECUTABLE})
  # the variables FLITE and ENGINE are set below to be used by the
  # CMake-based HTS voice building scripts
  if("${executable}" STREQUAL "flite")
    set(FLITE ${IBISENGINE_${executable}_EXECUTABLE})
  endif()
  if("${executable}" STREQUAL "hts_engine")
    set(ENGINE ${IBISENGINE_${executable}_EXECUTABLE})
  endif()
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  IbisEngine DEFAULT_MSG
  IBISENGINE_LIBRARY
  IBISENGINE_INCLUDE_DIR
  IBISENGINE_EXECUTABLE
  )

if(IBISENGINE_FOUND)
  set(IBISENGINE_INCLUDE_DIRS ${IBISENGINE_INCLUDE_DIR})
  set(IBISENGINE_LIBRARIES    ${IBISENGINE_LIBRARY})
  set(IBISENGINE_EXECUTABLES  ${IBISENGINE_EXECUTABLE})
endif()

mark_as_advanced(IBISENGINE_LIBRARY IBISENGINE_INCLUDE_DIR IBISENGINE_EXECUTABLE)
