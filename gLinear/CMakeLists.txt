# Copyright (c) 1995-2010, Stellenbosch University
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Stellenbosch University nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#
# gLinear CMakeLists.txt
#

message(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "Build  directory: $ENV{PWD}")

cmake_minimum_required(VERSION 3.13)

message(STATUS "Configuring gLinear ...")

# Allow -DCMAKE_BUILD_TYPE to be specified on the command line
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
else()
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
endif()

# declare project with default languages, i.e. C and C++
if("${PROJECT_TITLE}" STREQUAL "")
  set(PROJECT_TITLE gLinear)
endif()
project(${PROJECT_TITLE})

# define CMAKE_MODULE_PATH if it's not defined
if("${CMAKE_MODULE_PATH}" STREQUAL "")
  set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
endif()

# This ensures that TARGET_PLATFORM is defined
include(TargetPlatform)

# Build static and shared versions of libraries using CMake requires all the source
# files to be compiled twice as target names has to be unique.
option(BUILD_SHARED_LIBS "Set to OFF to build static gLinear" ON)
option(INDEX_CHECKING    "Set to ON to enable index checking (slow)" OFF)
# Choose whether or not to use boost.  At the time of writing, the only thing
# from boost we use are the pools.  Nothing there to link.
option(WANT_BOOST        "Set to ON to look for Boost support" OFF)
option(WANT_DOCS         "Set to ON to build PDF documentation" OFF)

set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-std=c++17 -O3 ") # -DTESTING_MODE -DDEBUG_MODE -DCHK_IDX -DDEBUG_MODE -pedantic
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ctor-dtor-privacy")

# C++ Compiler flags.
set(WARNINGS "-Wall -W -Wchar-subscripts -Wformat=2 -Wpointer-arith -Wwrite-strings -Wconversion -Wredundant-decls -Wctor-dtor-privacy -Wnon-virtual-dtor -Wcast-qual -Wold-style-cast -Woverloaded-virtual")

# Switch off this warning on Mac OS X, as its gcc seems confused by protected members
if(APPLE)
  set(WARNINGS "${WARNINGS} -Wno-ctor-dtor-privacy")
endif()

# -Weffc++ -Wshadow -Wfloat-equal

set(CMAKE_CXX_FLAGS                "-Wall -W -Wchar-subscripts -Wpointer-arith -Wcast-qual -Wwrite-strings -Wconversion -Wno-old-style-cast -Wctor-dtor-privacy -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS                "-std=c++17 ${CMAKE_CXX_FLAGS} -fPIC -ansi -pedantic -ftemplate-depth-99 ${WARNINGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++17 -O3 -g")
set(CMAKE_CXX_FLAGS_RELEASE        "-std=c++17 -O3 -g")
set(CMAKE_CXX_FLAGS_DEBUG          "-std=c++17 -g3 -ggdb -O0 -DDEBUG_MODE")

# Libraries required by libgLinear
set(gLinear_LIBS)
include(CheckLibraryExists)
include(CheckFunctionExists)

if(WANT_BOOST)
  find_package(Boost COMPONENTS system)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_definitions(-DHAVE_BOOST)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()

SET(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> rsc <TARGET> <OBJECTS>")

if(APPLE)
  if(DEFINED CMAKE_OSX_ARCHITECTURES)
    set(CMAKE_OSX_ARCHITECTURES ${CMAKE_OSX_ARCHITECTURES} CACHE STRING "Choose build architecture.")
  else()
    set(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "Choose build architecture.")
  endif()
endif()

set(gLinear_SRCS
  gError.cc
  gIndexRange.cc
  gLapackWrap.cc
  gLinearAlloc.cc
  gLinearMacros.cc
  matrix_funcs.cc
  )

set(gLinear_TSRCS
  gColVector.tcc gFortranTraits.tcc gLapackWrap.tcc gLinearAlgo.tcc
  gLinearFunctors.tcc gLinearMem.tcc gMatrix.tcc gMatrixBase.tcc
  gMatrixCore.tcc gMatrixIterator.tcc gMatrixMemBase.tcc
  gMatrixOpers.tcc gRowVector.tcc gVectorBase.tcc gVectorCore.tcc
  gVectorIterator.tcc gVectorMemBase.tcc gVectorOpers.tcc
  matrix_funcs.tcc vector_funcs.tcc
  )

set(gLinear_HSRCS
  Matrix.h Vector.h exceptions.h gColVector.h gCppTraits.h gError.h
  gFortranTraits.h gIndexRange.h gLapackHeaders.h gLapackWrap.h
  gLinear.h gLinearAlgo.h gLinearAlloc.h gLinearFunctors.h
  gLinearMacros.h gLinearMem.h gMatrix.h gMatrixBase.h gMatrixCore.h
  gMatrixIterator.h gMatrixMemBase.h gMatrixOpers.h gMatrixTraits.h
  gRowVector.h gVectorBase.h gVectorCore.h gVectorIterator.h
  gVectorMemBase.h gVectorOpers.h gVectorTraits.h matrix_funcs.h
  vector_funcs.h
  )

# Rule to build the gLinear library
set(gLinear_LIBNAME gLinear)

add_library(${gLinear_LIBNAME}
  ${gLinear_SRCS}
  ${gLinear_TSRCS}
  ${gLinear_HSRCS}
  )

# link to Accelerate.framework on OS X and liblapack elsewhere
if(APPLE)
  include(${CMAKE_SOURCE_DIR}/AddFramework.cmake)
  add_framework(${gLinear_LIBNAME} Accelerate)
else()
  check_library_exists(lapack zgesvd_ "" HAVE_LAPACK)
  if(HAVE_LAPACK)
    target_link_libraries(${gLinear_LIBNAME} lapack)
  else()
    message(FATAL_ERROR "can't find lapack")
  endif()
endif()

if(Boost_FOUND)
  target_link_libraries(${gLinear_LIBNAME} ${Boost_LIBRARIES} pthread)
endif()

# install library target
install(TARGETS ${gLinear_LIBNAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Rule to build the testnew executable
if(NOT CMAKE_CROSSCOMPILING)
  add_executable(testnew
    testnew.cc
    )
  target_link_libraries(testnew ${gLinear_LIBNAME})
endif()

# Documentation
if(WANT_DOCS)
  add_subdirectory(doc)
endif()
