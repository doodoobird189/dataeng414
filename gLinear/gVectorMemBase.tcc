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

/*
 * $Id$
 */

// standard headers
#include <cstddef>  // NULL
#include <new>  // placement new
#include <iostream>  // cout, endl


GLINEAR_NAMESPACE_BEGIN


template<typename T_Vector>
const char gVectorMemBase<T_Vector>::className[] = "gVectorMemBase<T_Vector>";



// no-init ctor
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(gLinearNoInit) :
  stor(NULL) {
}

// ctor specifying size and dflt ctor (if new_nElem == 0).
// alloc stor for elem elements. incr ref cnt.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(T_Index new_nElem) :
  nElem(new_nElem) {
  stor = new T_StorageElem(new_nElem);
  stor->addref();
  startAddr = stor->getStartAddr();
}

// copy ctor. alias arg gVectorMemBase vec. incr ref cnt.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(const gVectorMemBase<T_Vector>& vec) :
  nElem(vec.nElem),
  stor(vec.stor),
  startAddr(vec.startAddr) {
  stor->addref();
}


// ctor and init stor. incr ref cnt.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(T_Index new_nElem,
                                         const T_StorageElem* new_stor,
                                         T_ConstPointer new_startAddr) :
  nElem(new_nElem),
  stor( const_cast<T_StorageElem*>(new_stor) ),
  startAddr( const_cast<T_Pointer>(new_startAddr) ) {
  stor->addref();
}

// ctor from dense C/C++ array. Allocate elem elements and init from addr.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(T_ConstPointer addr,
                                         T_Index elem) :
  nElem(elem) {
  // ctruct and init
  stor = new T_StorageElem(addr, elem);
  stor->addref();
  startAddr = stor->getStartAddr();
}

// ctor from existing dense C/C++ array.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(T_ConstPointer addr,
                                         T_Index elem,
                                         const vector_no_copy_tag&) :
  nElem(elem) {
  // ctruct and init
  stor = new T_StorageElem();
  stor->setExisting(addr, elem);
  startAddr = const_cast<T_Pointer>( addr );
}

// ctor and init with arg.
template<typename T_Vector>
gVectorMemBase<T_Vector>::gVectorMemBase(T_ConstReference fill,
                                         T_Index elem) :
  nElem(elem) {
  // ctruct and init
  stor = new T_StorageElem(fill, elem);
  stor->addref();
  startAddr = stor->getStartAddr();
}

// dtor. decr ref cnt. deallocate RefCntStorage<T> object if ref cnt == 0.
template<typename T_Vector>
gVectorMemBase<T_Vector>::~gVectorMemBase(void) {
  if (stor) {
    stor->rmref();
  } // if
}

// configure
template<typename T_Vector>
void gVectorMemBase<T_Vector>::configure(T_Index new_nElem,
                                         const T_StorageElem* new_stor,
                                         T_ConstPointer new_startAddr) {
  // destruct and use the placement-new to call ctor to configure
  this->~gVectorMemBase();
  new(this) gVectorMemBase(new_nElem, new_stor, new_startAddr);
}


// dump info
template<typename T_Vector>
void gVectorMemBase<T_Vector>::specialInfo(int i) const {
  std::cout << "gVectorMemBase<T_Vector>" << std::endl;
  std::cout << "  nElem: " << nElem
       << " stor: " << stor
       << " startAddr: " << startAddr
       << std::endl;
  std::cout << " gVectorMemBase<T_Vector>::";
  stor->info(i);
}

// return nElem
template<typename T_Vector>
inline typename gVectorMemBase<T_Vector>::T_Index
gVectorMemBase<T_Vector>::specialSize(void) const {
  return nElem;
}

// resize
template<typename T_Vector>
inline bool gVectorMemBase<T_Vector>::specialResize(T_Index elem) {
  nElem = elem;
  stor->setSize(elem);
  startAddr = stor->getStartAddr();
  return true;
}

// return references
template<typename T_Vector>
inline typename gVectorMemBase<T_Vector>::T_Index
gVectorMemBase<T_Vector>::specialGetRefs(void) const {
  return stor->getRefs();
}

// deepen
template<typename T_Vector>
void gVectorMemBase<T_Vector>::specialDeepen(void) const {
  // save starting index.
  T_Index delta = static_cast<T_Index>(startAddr - stor->getStartAddr());
  // duplicate stor
  T_StorageElem* tmp = stor->duplicate();
  // decrement ref on old stor, reassign stor and inc ref on new stor
  stor->rmref();
  GLINEAR_CONST_OVERRIDE_THIS(gVectorMemBase<T_Vector>)->stor = tmp;
  stor->addref();
  // restore starting index
  GLINEAR_CONST_OVERRIDE_THIS(gVectorMemBase<T_Vector>)->startAddr = stor->getStartAddr() + delta;
}

// const alias
template<typename T_Vector>
inline void gVectorMemBase<T_Vector>::specialAlias(const gVectorMemBase<T_Vector>& vec) const {
  // set size
  GLINEAR_CONST_OVERRIDE_THIS(gVectorMemBase<T_Vector>)->nElem = vec.nElem;
  // decr ref on old stor, reassign stor and incr ref on new stor.
  // Note: this is not atomic and therefore not necessarily thread-safe.
  stor->rmref();
  GLINEAR_CONST_OVERRIDE_THIS(gVectorMemBase<T_Vector>)->stor = vec.stor;
  stor->addref();
  // set start address
  GLINEAR_CONST_OVERRIDE_THIS(gVectorMemBase<T_Vector>)->startAddr = vec.startAddr;
}

// alias
template<typename T_Vector>
inline void gVectorMemBase<T_Vector>::specialAlias(gVectorMemBase<T_Vector>& vec) {
  // set size
  nElem = vec.nElem;
  // decr ref on old stor, reassign stor and incr ref on new stor.
  // Note: this is not atomic and therefore not necessarily thread-safe.
  stor->rmref();
  stor = vec.stor;
  stor->addref();
  // set start address
  startAddr = vec.startAddr;
}

// return nElem
template<typename T_Vector>
inline typename gVectorMemBase<T_Vector>::T_Index
gVectorMemBase<T_Vector>::getNElem(void) const {
  return nElem;
}

// return ptr to stor
template<typename T_Vector>
inline typename gVectorMemBase<T_Vector>::T_StorageElem*
gVectorMemBase<T_Vector>::getStor(void) const {
  return stor;
}

// return startAddr
template<typename T_Vector>
inline typename gVectorMemBase<T_Vector>::T_Pointer
gVectorMemBase<T_Vector>::getStartAddr(void) const {
  return startAddr;
}

GLINEAR_NAMESPACE_END
