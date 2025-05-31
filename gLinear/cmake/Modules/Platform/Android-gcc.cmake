# Android CMake Compiler definition

# C Compile Flags
SET(CMAKE_SHARED_LIBRARY_C_FLAGS      "-fpic")     # -pic 
SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
SET(CMAKE_C_FLAGS_INIT                "-march=armv5te -mtune=xscale -msoft-float -mthumb-interwork -ffunction-sections -funwind-tables -fstack-protector -fno-short-enums -D__ARM_ARCH_5__ -D__ARM_ARCH_5T__ -D__ARM_ARCH_5E__ -D__ARM_ARCH_5TE__ -mthumb -Os -fomit-frame-pointer -fno-strict-aliasing -finline-limit=64 -O2 -DNDEBUG -g")
SET(CMAKE_C_FLAGS_DEBUG_INIT          "-fno-omit-frame-pointer -fno-strict-aliasing")
SET(CMAKE_C_FLAGS_MINSIZEREL_INIT     "-O2 -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300")
SET(CMAKE_C_FLAGS_RELEASE_INIT        "-O2 -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300")

# C++ Compile Flags
SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS      "-fpic")     # -pic 
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
SET(CMAKE_CXX_FLAGS_INIT                "${CMAKE_C_FLAGS_INIT} -fno-exceptions -fno-rtti")
SET(CMAKE_CXX_FLAGS_DEBUG_INIT          "-fno-omit-frame-pointer -fno-strict-aliasing")
SET(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "-fno-omit-frame-pointer -fno-strict-aliasing")
SET(CMAKE_CXX_FLAGS_RELEASE_INIT        "-O2 -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300")

# SET(CMAKE_SHARED_MODULE_C_FLAGS ${CMAKE_SHARED_LIBRARY_C_FLAGS})

# Misc flags
SET(CMAKE_INCLUDE_SYSTEM_FLAG_C   "-I")
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-I")

# C Standard Libraries
# TARGET_PREBUILT_SHARED_LIBRARIES from NDK
SET(CMAKE_C_STANDARD_LIBRARIES_INIT)   # reset
SET(CMAKE_CXX_STANDARD_LIBRARIES_INIT) # reset
# if compiling static libraries use static standard libraries
IF(BUILD_SHARED_LIBS)
  SET(ANDROID_STDLIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
ELSE(BUILD_SHARED_LIBS)
  SET(ANDROID_STDLIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDIF(BUILD_SHARED_LIBS)

FOREACH(library c stdc++ m) # libthread_db.so
 SET(CMAKE_C_STANDARD_LIBRARIES_INIT   "${CMAKE_C_STANDARD_LIBRARIES_INIT} ${SYSROOT}/usr/lib/lib${library}${ANDROID_STDLIB_SUFFIX}")
 # SET(CMAKE_C_STANDARD_LIBRARIES_INIT   "${CMAKE_C_STANDARD_LIBRARIES_INIT} -l${library}")
ENDFOREACH(library)
# C++ Standard Libraries
FOREACH(library c stdc++ m) # libthread_db.so
  SET(CMAKE_CXX_STANDARD_LIBRARIES_INIT "${CMAKE_CXX_STANDARD_LIBRARIES_INIT} ${SYSROOT}/usr/lib/lib${library}${ANDROID_STDLIB_SUFFIX}")
  # SET(CMAKE_CXX_STANDARD_LIBRARIES_INIT "${CMAKE_CXX_STANDARD_LIBRARIES_INIT} -l${library}")
ENDFOREACH(library)

# Android's C Runtime objects
# NDK's TARGET_CRTBEGIN_STATIC_O TARGET_CRTBEGIN_DYNAMIC_O TARGET_CRTEND_O
# NDK-1.6 does not use crtbegin_static.o so it's not included here
SET(TARGET_CRTBEGIN_DYNAMIC_O "${SYSROOT}/usr/lib/crtbegin_dynamic.o")
SET(TARGET_CRTEND_O           "${SYSROOT}/usr/lib/crtend_android.o")
#SET(TARGET_CRTBEGIN_DYNAMIC_O "")
# SET(TARGET_CRTEND_O           "")

execute_process(
  COMMAND ${CMAKE_C_COMPILER} -mthumb-interwork -print-libgcc-file-name
  OUTPUT_VARIABLE ANDROID_LIBGCC
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# to prevent reloc_library[1172]:  cannot locate '__aeabi_fadd'... error
SET(ANDROID_LIBGCC "-Wl,--no-undefined ${CMAKE_SHARED_LIBRARY_RPATH_LINK_C_FLAG}${SYSROOT}/usr/lib ${ANDROID_LIBGCC}") # -Wl,--whole-archive

SET(CMAKE_C_CREATE_SHARED_LIBRARY
  "<CMAKE_C_COMPILER> -nostdlib <CMAKE_SHARED_LIBRARY_SONAME_C_FLAG><TARGET_SONAME> <LINK_FLAGS> <OBJECTS> <LINK_LIBRARIES> ${ANDROID_LIBGCC} -o <TARGET>")

SET(CMAKE_C_CREATE_SHARED_MODULE ${CMAKE_C_CREATE_SHARED_LIBRARY})

SET(CMAKE_C_CREATE_STATIC_LIBRARY
  "<CMAKE_AR> crs <TARGET> <OBJECTS>")

SET(CMAKE_C_COMPILE_OBJECT
  "<CMAKE_C_COMPILER> <DEFINES> <FLAGS> -c <SOURCE> -o <OBJECT>")

SET(CMAKE_C_LINK_EXECUTABLE
  "<CMAKE_C_COMPILER> -nostdlib <LINK_FLAGS> ${TARGET_CRTBEGIN_DYNAMIC_O} <OBJECTS> <LINK_LIBRARIES> ${ANDROID_LIBGCC} ${TARGET_CRTEND_O} -o <TARGET>")

# reuse C rules for C++
SET(CMAKE_CXX_CREATE_SHARED_LIBRARY "${CMAKE_C_CREATE_SHARED_LIBRARY}")
SET(CMAKE_CXX_CREATE_SHARED_MODULE  "${CMAKE_C_CREATE_SHARED_MODULE}")
SET(CMAKE_CXX_CREATE_STATIC_LIBRARY "${CMAKE_C_CREATE_STATIC_LIBRARY}")
SET(CMAKE_CXX_COMPILE_OBJECT        "${CMAKE_C_COMPILE_OBJECT}")
SET(CMAKE_CXX_LINK_EXECUTABLE       "${CMAKE_C_LINK_EXECUTABLE}")
