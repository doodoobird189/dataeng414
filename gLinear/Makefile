
# Copyright (c) 1995-2006, Stellenbosch University
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
# gLinear Makefile
#
# Targets:
#
# all:	   Build library and test application.
# lib:     Build library.
# testnew: Build test application.
# docs:    Generate Doxygen docs.
# clean:   Remove object, core, executable and library files.
#

# Set this to 1 to tweak the Makefile for Mac OS X
macosx=0

# Choose whether or not to use boost.  At the time of writing, the only thing
# from boost we use are the pools.  Nothing there to link.
use_boost=0

# C++ compiler executable
CXX=g++ # -3.4

# debug build? else you get a mean and lean optimised build
debug=0

# C++ Compiler flags.
WARNINGS=-Wall -W -Wchar-subscripts -Wformat=2 -Wpointer-arith \
         -Wwrite-strings -Wconversion -Wredundant-decls \
         -Wctor-dtor-privacy -Wnon-virtual-dtor -Wcast-qual \
         -Wold-style-cast -Woverloaded-virtual
# Switch off this warning on Mac OS X, as its gcc seems confused by protected members
ifeq (1,${macosx})
WARNINGS += -Wno-ctor-dtor-privacy
endif
# -Weffc++ -Wshadow -Wfloat-equal
CXXFLAGS=-fPIC -ansi -pedantic -ftemplate-depth-99 $(WARNINGS)

ifeq (1,${debug})
# Uncomment the defines below to get even more debugging aid. this
# dramatically slows down the code.
CXXFLAGS += -ggdb3 # -DDEBUG_MODE # -DCHK_IDX
OBJSFX = .od
SHAREDLIBSFX = d.so
STATICLIB = d.a
EXESFX = d
else
CXXFLAGS += -O3
OBJSFX = .o
SHAREDLIBSFX = .so
STATICLIB = .a
EXESFX = 
endif

ifeq (1,${macosx})
SHAREDLIBSFX=$(SHAREDLIBSFX:.so=.dylib)
endif

BOOST_ROOT = $(HOME)/opt
BOOST_INCLUDE = $(BOOST_ROOT)/include/boost-1_33_1
BOOST_LIB = $(BOOST_ROOT)/lib

ifeq (1,${use_boost})
CXXFLAGS += -I$(BOOST_INCLUDE) -DHAVE_BOOST -pthread
endif

# Libraries required by libgLinear
ifeq (1,${macosx})
LIBS=-framework veclib
else
LIBS=-llapack -lblas -lgfortran
endif

ifeq (1,${use_boost})
LIBS+=-lpthread
endif

# Directory where Doxygen HTML documentation is created
DOCDIR=./html
# Doxygen command to use
DOCGEN=/usr/bin/doc++  -v  -p  -f -d $(DOCDIR)

# List of files to remove when executing "clean" target.
RMLIST= $(CPPSOURCES:.cc=.o) \
	$(CPPSOURCES:.cc=.od) \
	$(CPPSOURCES:.cc=.dep) \
	libgLinear*.a \
	libgLinear*.so \
	testnew.o \
	testnew \
	test.out \
	core

# List of source files
SRC=gLinear.h \
	gError.h gError.cc \
	gIndexRange.h gIndexRange.cc \
	gLinearMacros.h gLinearMem.h \
	gVectorTraits.h \
	gVectorBase.h gVectorBase.tcc \
	gVectorMemBase.h gVectorMemBase.tcc \
	gVectorCore.h \
	gVectorIterator.h gVectorIterator.tcc \
	gRowVector.h gColVector.h \
	gVectorOpers.h gVectorOpers.tcc \
	gMatrixTraits.h \
	gMatrixBase.h gMatrixBase.tcc \
	gMatrixMemBase.h gMatrixMemBase.tcc \
	gMatrixCore.h \
	gMatrixIterator.h gMatrixIterator.tcc \
	gMatrix.h \
	gMatrixOpers.h gMatrixOpers.tcc \
	gLinearFunctors.h \
	gLapackWrap.h gLapackWrap.tcc gLapackWrap.cc gLinearAlloc.cc \
	vector_funcs.h vector_funcs.tcc \
	matrix_funcs.h matrix_funcs.cc matrix_funcs.tcc

CPPSOURCES=gError.cc \
	gIndexRange.cc \
	gLapackWrap.cc \
	gLinearAlloc.cc \
	gLinearMacros.cc \
	matrix_funcs.cc


all:   lib testnew$(EXESFX)

testnew$(EXESFX): testnew$(OBJSFX) lib
	$(CXX) -o $@ $< libgLinear$(STATICLIB) $(LIBS)

ifeq (1,${macosx})
lib:	libgLinear$(STATICLIB) libgLinear.dylib
else
lib:	libgLinear$(STATICLIB) libgLinear$(SHAREDLIBSFX)
endif

libgLinear$(STATICLIB): $(CPPSOURCES:.cc=$(OBJSFX))
	ar rscv $@ $^

# Linux shared library
libgLinear$(SHAREDLIBSFX): $(CPPSOURCES:.cc=$(OBJSFX))
	$(CXX) -shared -o $@ $^ $(LIBS)

# Mac OS X shared library
libgLinear.dylib: $(CPPSOURCES:.cc=$(OBJSFX))
	$(CXX) -dynamiclib -o $@ $^ $(LIBS)


#
# Include all dependency lists here, unless we are 'clean'ing or making 'docs'
#
# This will ensure the generation of .dep files from source files
#
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),docs)
# -include $(CSOURCES:.c=.dep)
-include $(CPPSOURCES:.cc=.dep)
endif
endif

#
# These rules generate the automatic dependency lists, using the compiler
#
# This follows the example in the GNU make manual, section 4.12
#
%.dep: %.c
	@ echo generating $@ . . .
	@ $(SHELL) -ec '$(CC) -MM $(CFLAGS) $< \
                        | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
                        [ -s $@ ] || rm -f $@'

%.dep: %.cc
	@ echo generating $@ ...
	@ $(SHELL) -ec '$(CC) -MM $(CFLAGS) $< \
                        | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
                        [ -s $@ ] || rm -f $@'


%.o: %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

%.od: %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean:
	for file in $(RMLIST); do \
		[ -f $$file ] && rm $$file; \
	done; \
        rm -f *~; \
	rm -rf $(DOCDIR); \

docs:
	[ -d $(DOCDIR) ] && rm -rf $(DOCDIR)
	$(DOCGEN) docroot.dxx;

.PHONY:	docs
