# select c++ compiler
set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-std=c++17 -O3 ") # -DDEBUG_MODE -DCHK_IDX -DDEBUG_MODE -pedantic

# to enable warnings
add_definitions ("-Wall")

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#
# The OTHER .cc source files
set(PRLITE_OBJECTS
  prlite_alloc.cc prlite_colourmap.cc prlite_endian.cc
  prlite_fnamepar.cc prlite_genmat.cc prlite_logging.cc
  prlite_matrix.cc prlite_mrandom.cc
  prlite_my_any.cc prlite_shuffles.cc prlite_stdfun.cc
  prlite_stringlist.cc prlite_stringops.cc prlite_testing.cc
  prlite_thread.cc prlite_zfstream3.cc
  )

# prlite_mean_cov.cc

#
# The OTHER .hpp header files
set(PRLITE_HEADERS
  prlite_alloc.hpp prlite_colourmap.hpp prlite_endian.hpp
  prlite_error.hpp prlite_fnamepar.hpp prlite_gendefs.hpp
  prlite_genmat.hpp prlite_genvec.hpp prlite_logging.hpp
  prlite_matflips.hpp prlite_matrix.hpp prlite_mrandom.hpp
  prlite_my_any.hpp prlite_outerquad_l2u.hpp prlite_rcptr.hpp
  prlite_sorting.hpp prlite_stdfun.hpp prlite_shuffles.hpp
  prlite_stlvecs.hpp prlite_stringlist.hpp prlite_stringops.hpp
  prlite_testing.hpp prlite_thread.hpp prlite_vector.hpp
  prlite_zfstream.hpp prlite_zfstream3.hpp
  )
# prlite_mean_cov.hpp

# The .tcc template source files
set(PRLITE_TEMPLATES
    prlite_endian.tcc
    prlite_genmat.tcc prlite_genvec.tcc prlite_matflips.tcc prlite_matrix.tcc
    prlite_outerquad_l2u.tcc prlite_rcptr.tcc prlite_sorting.tcc prlite_stdfun.tcc
    prlite_stlvecs.tcc prlite_vector.tcc
  )
# prlite_mean_cov.tcc
#
# Create a variable called PRLITE_SOURCES containing all code files:
set(PRLITE_SOURCES ${PRLITE_OBJECTS} ${PRLITE_HEADERS} ${PRLITE_TEMPLATES})
#

include_directories(${GLINEAR_INCLUDE_DIRS} )

# this builds the prlite library (libprlite.so)
add_library(prlite SHARED ${PRLITE_SOURCES})
target_link_libraries(prlite ${GLINEAR_LIBRARIES}  ssl crypto z)
