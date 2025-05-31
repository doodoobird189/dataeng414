
set(HTS_POSSIBLE_BIN_PATHS
  /usr/local/HTS-2.3alpha/bin
  /usr/local/HTS-2.2/bin
  /usr/local/HTS-2.2beta/bin
  ${HTS_SEARCH_PATH}
  )

INCLUDE(FindPackageHandleStandardArgs)

if(HTS_FIND_COMPONENTS)
  foreach(component ${HTS_FIND_COMPONENTS})
    string(TOUPPER ${component} _COMPONENT)
    # message(STATUS "HTS component: ${_COMPONENT}")
    find_program(${_COMPONENT}
      NAMES ${component}
      PATHS ${HTS_POSSIBLE_BIN_PATHS}
      )
    if(${_COMPONENT})
      FIND_PACKAGE_HANDLE_STANDARD_ARGS(HTS DEFAULT_MSG ${_COMPONENT})
    else(${_COMPONENT})
      message(FATAL_ERROR "HTS: can't find ${component}")
    endif(${_COMPONENT})
  endforeach(component ${HTS_FIND_COMPONENTS})
endif(HTS_FIND_COMPONENTS)



