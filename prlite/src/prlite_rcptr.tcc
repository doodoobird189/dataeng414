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
  $Id$
*/


// standard headers
#include <sstream>  // ostringstream
#include <iostream>  // endl, ostream
#ifdef ABORT_NOT_THROW
#include <cstdlib>  // abort
#endif // ABORT_NOT_THROW

// patrec headers
#include "prlite_alloc.hpp"

 //namespace prlite{

template<typename T>
inline RCPtr<T>::RCPtr()
    : oPtr(0) {
  cPtr = refCountNew();
  *cPtr = 1;
} // default constructor

template<typename T>
inline RCPtr<T>::RCPtr(const RCPtr& rcp)
    : oPtr(rcp.oPtr),
      cPtr(rcp.cPtr) {
  (*cPtr)++;
} // copy constructor

template<typename T>
template<typename FromType>
inline RCPtr<T>::RCPtr(const RCPtr<FromType>& rcp)
    : oPtr( rcp.objPtr() ),
      cPtr( rcp.cntPtr() ) {
  (*cPtr)++;
} // constructor

template<typename T>
inline RCPtr<T>::~RCPtr() {
  if ( !( --(*cPtr) ) ) {
    delete oPtr;
    oPtr = 0;
    refCountDelete(cPtr);
    cPtr = 0;
  } // if
} // destructor

template<typename T>
inline RCPtr<T>::RCPtr(T* justNewed)
    : oPtr(justNewed) {
  cPtr = refCountNew();
  *cPtr = 1;
} // constructor

template<typename T>
inline RCPtr<T>::RCPtr(T* objectPtr,
                       unsigned* countPtr)
    : oPtr(objectPtr),
      cPtr(countPtr) {
  (*cPtr)++;
} // constructor

template<typename T>
inline RCPtr<T>& RCPtr<T>::operator=(const RCPtr& rhs) {
  return assignFrom( rhs.objPtr(), rhs.cntPtr() );
} // operator=

template<typename T>
template<typename FromType>
inline RCPtr<T>& RCPtr<T>::operator=(const RCPtr<FromType>& rhs) {
  return assignFrom( rhs.objPtr(), rhs.cntPtr() );
} // operator=

template<typename T>
inline T* RCPtr<T>::operator->() const {

  // adt //
  // 05/04/2001 //
  // 17:16 //
  // Added a check for dereference of a null object pointer

#ifdef DEBUG_MODE
  if (!oPtr) {
    #ifdef ABORT_NOT_THROW
    std::cerr << '\n' << "ERROR => " << "RCPtr::operator-> : Dereferencing null pointer\n" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    std::abort();
    #else
    std::ostringstream errMsg;
    errMsg << '\n' << "ERROR => " << "RCPtr::operator-> : Dereferencing null pointer\n" << std::endl;
    errMsg << __FILE__ << ":" << __LINE__ << std::endl;
    throw errMsg.str();
    #endif // ABORT_NOT_THROW
  } // if
#endif // DEBUG_MODE

  return oPtr;

} // operator->

template<typename T>
inline T& RCPtr<T>::operator*() const {

  // adt //
  // 05/04/2001 //
  // 17:16 //
  // Added a check for dereferencing a null object pointer

#ifdef DEBUG_MODE
  if (!oPtr) {
    #ifdef ABORT_NOT_THROW
    std::cerr << '\n' << "ERROR => " << "RCPtr::operator-> : Dereferencing null pointer\n" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    std::abort();
    #else
    std::ostringstream errMsg;
    errMsg << '\n' << "ERROR => " << "RCPtr::operator-> : Dereferencing null pointer\n" << std::endl;
    errMsg << __FILE__ << ":" << __LINE__ << std::endl;
    throw errMsg.str();
    #endif // ABORT_NOT_THROW
  } // if
#endif // DEBUG_MODE

  return *oPtr;

} // operator*

template<typename T>
inline bool RCPtr<T>::operator<(const RCPtr& rhs) const {
  return (oPtr < rhs.oPtr);
} // operator<

template<typename T>
inline bool RCPtr<T>::operator==(const RCPtr& rhs) const {
  return (oPtr == rhs.oPtr);
} // operator<

template<typename T>
inline bool RCPtr<T>::operator!=(const RCPtr& rhs) const {
  return (oPtr != rhs.oPtr);
} // operator<

template<typename T>
inline RCPtr<T>::operator bool() const {
  return oPtr != 0;
} // operator bool

//=====================others======================================

template<typename T>
inline T* RCPtr<T>::objPtr() const {
  return oPtr;
} // objPtr

template<typename T>
inline unsigned* RCPtr<T>::cntPtr() const {
  return cPtr;
} // cntPtr

template<typename T>
inline void RCPtr<T>::makeImmortal() {
  (*cPtr)++;
} // makeImmortal

template<typename T>
inline bool RCPtr<T>::shared() const {
  return oPtr && *cPtr > 1;
} // shared

template<typename T>
inline RCPtr<T>& RCPtr<T>::assignFrom(T* objectPtr,
                                      unsigned* countPtr) {
  (*cPtr)--;
  (*countPtr)++;
  if ( !(*cPtr) ) {
    delete oPtr;
    refCountDelete(cPtr);
  } // if
  cPtr = countPtr;
  oPtr = objectPtr;
  return *this;
} // assignFrom

//=====================friends=====================================

template<typename ToType, typename FromType>
RCPtr<ToType> dynamicCast(const RCPtr<FromType>& fromPtr) {
  ToType* rawPtr = dynamic_cast<ToType*>( fromPtr.objPtr() );
  if (rawPtr) {
    return RCPtr<ToType>( rawPtr, fromPtr.cntPtr() );
  } // if
  else {
    return RCPtr<ToType>(0);
  } // else
} // dynamicCast

template<typename ToType, typename FromType>
RCPtr<ToType> constCast(const RCPtr<FromType>& fromPtr) {
  ToType* rawPtr = const_cast<ToType*>( fromPtr.objPtr() );
  return RCPtr<ToType>( rawPtr, fromPtr.cntPtr() );
} // constCast

template<typename ToType, typename FromType>
bool assignTo(RCPtr<ToType>& toPtr,
              const RCPtr<FromType>& fromPtr) {
  ToType* rawPtr = dynamic_cast<ToType*>( fromPtr.objPtr() );
  if (rawPtr) {
    return toPtr.assignFrom( rawPtr, fromPtr.cntPtr() );
  } // if
  else {
    return toPtr = RCPtr<ToType>(0);
  } // else
} // assignTo

template<typename T>
inline std::ostream& operator<<(std::ostream& file,
        const RCPtr<T>& rcp) {
  return file << (*rcp.cPtr) << " refs to " << rcp.oPtr;
} // operator<<

//} // namespace prlite
