# Android Platform CMake definition

# GCC is the default compiler on Android
# SET(CMAKE_DL_LIBS "dl")
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG     "-Wl,-rpath,")
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP ":")
SET(CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG  "-Wl,-rpath-link=")
SET(CMAKE_SHARED_LIBRARY_SONAME_C_FLAG      "-Wl,-soname,")
SET(CMAKE_SHARED_LIBRARY_SONAME_CXX_FLAG    "-Wl,-soname,")

# Shared libraries with no builtin soname may not be linked safely by
# specifying the file path.
SET(CMAKE_PLATFORM_USES_PATH_WHEN_NO_SONAME 1)

# Disable rpath
# SET(CMAKE_SKIP_RPATH "YES" CACHE BOOL)
SET(CMAKE_LIBRARY_PATH_FLAG "-L")
# SET(CMAKE_LIBRARY_PATH_FLAG "${CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG}${SYSROOT/usr/lib}")

# Executable Linker Flags
SET(CMAKE_EXE_LINKER_FLAGS_INIT "-Bdynamic -Wl,-dynamic-linker,/system/bin/linker -Wl,--gc-sections -Wl,-z,nocopyreloc")

# Shared Library Linker Flags
SET(CMAKE_SHARED_LINKER_FLAGS_INIT "-Wl,-shared,-Bsymbolic -Wl,--no-whole-archive -Wl,--error-unresolved-symbols -Wl,-z,defs")

# Module Linker Flags, same as for Shared Library
SET(CMAKE_MODULE_LINKER_FLAGS_INIT   "${CMAKE_SHARED_LINKER_FLAGS_INIT}")

# INCLUDE(Platform/UnixPaths)

# include the gcc (arm-eabi) compiler information
INCLUDE(Platform/Android-gcc)

SET(ANDROID 1)
SET(ANDROID_VERSION ${NDKVER})
