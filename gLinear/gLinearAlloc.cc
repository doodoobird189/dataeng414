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

// define HAVE_BOOST and point your compiler to the boost includes to
// use a special pool for allocating these.  The boost::pool library
// is a headers only library.


// in debug mode we don't want the boost pools, because the hide memory
// bugs.
#ifdef DEBUG_MODE
#undef HAVE_BOOST
#endif

#include "gLinearAlloc.h"
#include "gLinearMem.h"
// #include <iostream>
#include <exception>

#ifdef HAVE_BOOST
#include <boost/pool/singleton_pool.hpp>
#endif


GLINEAR_NAMESPACE_BEGIN


#ifdef HAVE_BOOST

struct __RefCntStorage_float_pool_tag { } ;
typedef boost::singleton_pool<
  __RefCntStorage_float_pool_tag,
  sizeof(RefCntStorage<double>),
  boost::default_user_allocator_new_delete,
  boost::details::pool::default_mutex,
  10*1024 // make it big, we use many.
  > __RefCntStorage_float_pool;

struct __RefCntStorage_double_pool_tag { } ;
typedef boost::singleton_pool<
  __RefCntStorage_double_pool_tag,
  sizeof(RefCntStorage<double>),
  boost::default_user_allocator_new_delete,
  boost::details::pool::default_mutex,
  10*1024 // make it big, we use many.
  > __RefCntStorage_double_pool;

#endif // HAVE_BOOST


void* allocRefCntStorageDouble() {
#ifdef HAVE_BOOST
  // std::cout << "allocRefCntStorageDouble() - pool" << std::endl;
  void* ret = __RefCntStorage_double_pool::malloc();
  if (!ret)
    throw std::bad_alloc();
  else
    return ret;
#else
  // std::cout << "allocRefCntStorageDouble() - new" << std::endl;
  return ::operator new( sizeof(RefCntStorage<double>) );
#endif
}

void freeRefCntStorageDouble(void* stor) {
#ifdef HAVE_BOOST
  // std::cout << "freeRefCntStorageDouble() - pool" << std::endl;
  __RefCntStorage_double_pool::free( stor );
#else
  // std::cout << "freeRefCntStorageDouble() - delete" << std::endl;
  ::operator delete ( reinterpret_cast<void*>(stor) );
#endif
}

void* allocRefCntStorageFloat() {
#ifdef HAVE_BOOST
  // std::cout << "allocRefCntStorageFloat() - pool" << std::endl;
  void* ret = __RefCntStorage_float_pool::malloc();
  if (!ret)
    throw std::bad_alloc();
  else
    return ret;
#else
  // std::cout << "allocRefCntStorageFloat() - new" << std::endl;
  return ::operator new( sizeof(RefCntStorage<float>) );
#endif
}

void freeRefCntStorageFloat(void* stor) {
#ifdef HAVE_BOOST
  // std::cout << "freeRefCntStorageFloat() - pool" << std::endl;
  __RefCntStorage_float_pool::free( stor );
#else
  // std::cout << "freeRefCntStorageFloat() - delete" << std::endl;
  ::operator delete ( reinterpret_cast<void*>(stor) );
#endif
}

GLINEAR_NAMESPACE_END


