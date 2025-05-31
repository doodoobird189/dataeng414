# Try to find PortAudio
# Once done, this will define
#
# PORTAUDIO_FOUND - system has PortAudio
# PORTAUDIO_INCLUDE_DIRS - the PortAudio include directories
# PORTAUDIO_LIBRARIES - link these to use PortAudio

if(PORTAUDIO_INCLUDE_DIRS AND PORTAUDIO_LIBRARIES)
    set(PORTAUDIO_FIND_QUIETLY TRUE)
endif()

# include dir
find_path(
  PORTAUDIO_INCLUDE_DIR portaudio.h
  PATHS /usr/local/include /usr/include
  )

# finally the library itself
find_library(
  PORTAUDIO_LIBRARY NAMES portaudio
  PATHS /usr/local/lib /usr/lib
  )

# Seems to run OK without libportaudiocpp so do not pull additional dependency
# find_library(libPortAudioCpp NAMES portaudiocpp)
# set(PORTAUDIO_LIBRARIES ${libPortAudio} ${libPortAudioCpp})
# handle the QUIETLY and REQUIRED arguments and set PORTAUDIO_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  PortAudio DEFAULT_MSG
  PORTAUDIO_LIBRARY PORTAUDIO_INCLUDE_DIR
  )

set(PORTAUDIO_LIBRARIES ${PORTAUDIO_LIBRARY})
set(PORTAUDIO_INCLUDE_DIRS ${PORTAUDIO_INCLUDE_DIR})

mark_as_advanced(PORTAUDIO_LIBRARY PORTAUDIO_INCLUDE_DIR)
