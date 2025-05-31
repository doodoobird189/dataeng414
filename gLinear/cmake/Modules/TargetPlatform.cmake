# TARGET_PLATFORM is always specified by the ./cmake/configure script
# If cmake is used without the configure script we need to make sure
# TARGET_PLATFORM is still defined. The value of TARGET_PLATFORM is
# used in the ./cmake/Modules/Find*.cmake scripts to find built
# dependencies, so make sure you adopt the same build directory
# structure as that produced by the configure script, unless of course
# you know what you're doing. For cross-compilation it is recommended
# to always use the configure script.
if("${TARGET_PLATFORM}" STREQUAL "")
  execute_process(COMMAND uname -s OUTPUT_VARIABLE TARGET_OS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND uname -m OUTPUT_VARIABLE TARGET_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(TARGET_PLATFORM "${TARGET_OS}-${CMAKE_BUILD_TYPE}-${TARGET_ARCH}")
endif()
message(STATUS "Target platform: ${TARGET_PLATFORM}")
