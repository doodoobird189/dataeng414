/*
Copyright (c) 1995-2006, Stellenbosch University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Stellenbosch University nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id$

// Author: A.S. Visagie

// Purpose: special functions for allocating the two most common types
// of RefCntStorage<>.  The reasoning for factoring their allocation
// out here is
// 1. so that we can control its instantiation: one DLL, and
// 2. we can replace the allocation with a pool and keep the fact
// complete hidden from all code except gLinearAlloc.cc

#ifndef GLINEARALLOC_HPP
#define GLINEARALLOC_HPP

#include <new>  // bad_alloc
#include "gLinearMacros.h"

GLINEAR_NAMESPACE_BEGIN

template<typename T>
class RefCntStorage;

GLINEAR_API void* allocRefCntStorageDouble();
GLINEAR_API void freeRefCntStorageDouble(void* stor);

GLINEAR_API void* allocRefCntStorageFloat();
GLINEAR_API void freeRefCntStorageFloat(void* stor);

GLINEAR_NAMESPACE_END

#endif // GLINEARALLOC_HPP
