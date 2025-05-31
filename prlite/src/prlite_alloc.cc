char ALLOC_CC_VERSION[] = "$Id$";

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

// in debug mode we don't want the boost pools, because the hide memory
// bugs.
#ifdef DEBUG_MODE
#undef HAVE_BOOST
#endif

#ifdef HAVE_BOOST
#include <boost/pool/singleton_pool.hpp>

struct __rcptr_refcount_pool_tag { } ;
typedef boost::singleton_pool<
  __rcptr_refcount_pool_tag,
  sizeof(unsigned),
  boost::default_user_allocator_new_delete,
  boost::details::pool::default_mutex,
  512/sizeof(unsigned) // make the NextSize slightly bigger than the
           // default, but still small enough for the
           // small objects heap.
  > __rcptr_refcount_pool;
#endif

#include "prlite_alloc.hpp"
#include <exception>

/** allocate a new unsigned int on the heap.  If BOOST is available,
    use a singleton pool.  This has a small speed advantage, and a
    huge memory fragmentation advantage. */
unsigned* refCountNew() {
#ifdef HAVE_BOOST
  unsigned* ret = reinterpret_cast<unsigned*>( __rcptr_refcount_pool::malloc() );
  if (!ret)
    throw std::bad_alloc();
  else
    return ret;
#else
  return new unsigned;
#endif
}

/** deallocate the unsigned int. */
void refCountDelete(unsigned* cPtr) {
#ifdef HAVE_BOOST
    __rcptr_refcount_pool::free( reinterpret_cast<void*>(cPtr) );
#else
    delete cPtr;
#endif
}
