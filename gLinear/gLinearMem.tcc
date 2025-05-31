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

// glinear headers
#include "exceptions.h"

// standard headers
#include <cstddef>  // NULL
#include <new>  // bad_alloc
#include <iostream>  // cout, endl
#include <sstream>  // ostringstream
#include <limits>  // numeric_limits
#ifdef ABORT_NOT_THROW
#include <cstdlib>  // abort
#endif // ABORT_NOT_THROW

GLINEAR_NAMESPACE_BEGIN

template<typename T>
const char RefCntStorage<T>::ClassName[] = "RefCntStorage<T>";


template<typename T>
inline RefCntStorage<T>::RefCntStorage(void) :
  refs(0),
  size(0),
  startAddr(NULL) {
} // default constructor

template<typename T>
inline RefCntStorage<T>::RefCntStorage(T_Index new_size) :
  refs(0),
  size(new_size) {
  try {
    startAddr = new T_Value[new_size];
  } catch (const std::bad_alloc&) {
    gErrorBadAlloc("RefCntStorage",
                   "RefCntStorage(T_Index new_size)",
                   "new() failed",
                   sizeof(T_Value),
                   new_size);
  } // catch
} // constructor

template<typename T>
inline RefCntStorage<T>::RefCntStorage(const T_Pointer addr,
               T_Index new_size) :
  refs(0),
  size(new_size) {
  startAddr = new T_Value[new_size];
  for (T_Index indx = 0; indx < new_size; indx++) {
    startAddr[indx] = addr[indx];
  } // for
} // constructor

template<typename T>
inline RefCntStorage<T>::RefCntStorage(const T_Value& fill,
               T_Index new_size) :
  refs(0),
  size(new_size) {
  startAddr = new T[new_size];
  T_Index indx = new_size;
  while (--indx) {
    startAddr[indx] = fill;
  } // while
} // constructor

template<typename T>
inline RefCntStorage<T>::~RefCntStorage(void) {
  if (startAddr != NULL &&
      refs != std::numeric_limits<T_Index>::max()) {
    delete[] startAddr;
  } // if
} // destructor

template<typename T>
inline void RefCntStorage<T>::info(int i) const {
  std::cout << "RefCntStorage<T>" << std::endl;
  std::cout << "  refs: " << refs
       << " size: " << size
       << " startAddr: " << startAddr
       << std::endl;
} // info

template<typename T>
inline bool RefCntStorage<T>::rmref(void) {

  if (refs == std::numeric_limits<T_Index>::max())
    return false; // do nothing

  switch (refs) {
    // pretty serious. This should not ever happen.
  case 0: {
    std::ostringstream errMsg;
    errMsg << '\n' << "ERROR => " << "RefCntStorage<T_Value>::rmref(void) -- Internal integrity failure detected." << std::endl;
    errMsg << __FILE__ << ":" << __LINE__ << std::endl;
#ifdef ABORT_NOT_THROW
    std::cerr << errMsg.str() << std::flush;
    std::abort();
#else
    throw GLinearError( errMsg.str() );
#endif // ABORT_NOT_THROW
    break;
  } // case 0:

  case 1:
    refs--;
    delete this;
    return true;
    break;
  default:
    refs--;
    break;
  } // switch
  return false;
} // rmref

template<typename T>
inline void RefCntStorage<T>::addref(void) {
  if (refs == std::numeric_limits<T_Index>::max())
    return; // do nothing
  else
    refs++;
} // addref

template<typename T>
inline T_Index RefCntStorage<T>::getSize(void) {
  return size;
} // getSize

template<typename T>
inline void RefCntStorage<T>::setSize(T_Index new_size) {
  // if this memory belongs to someone else, you cannot resize it.
  if (refs == std::numeric_limits<T_Index>::max()) {
    gErrorAliases("RefCntStorage",
      "setSize(T_Index new_size)",
      "Cannot resize vector when the memory does not belong to me.");
  }
  // free old storage
  if (startAddr != NULL) {
    delete[] startAddr;
  } // if
  // allocate new
  startAddr = new T_Value[new_size];
  size = new_size;
} // setSize

template<typename T>
inline T_Index RefCntStorage<T>::getRefs(void) {
  return refs;
} // getRefs

template<typename T>
inline T* RefCntStorage<T>::getStartAddr(void) {
  return startAddr;
} // getStartAddr

template<typename T>
inline RefCntStorage<T>* RefCntStorage<T>::duplicate(void) {
  RefCntStorage<T_Value>* tmp = new RefCntStorage<T_Value>(size);
  for (T_Index indx = 0; indx < size; indx++) {
    (tmp->startAddr)[indx] = startAddr[indx];
  } // for
  return tmp;
} // duplicate


// generic one, uses new
template<typename T>
inline void* RefCntStorage<T>::operator new(size_t size) {
  return ::operator new( size );
}

// generic one, uses delete
template<typename T>
inline void RefCntStorage<T>::operator delete(void* stor) {
  ::operator delete( stor );
}

// specialisation for double
template<>
inline void* RefCntStorage<double>::operator new(size_t) {
  return allocRefCntStorageDouble();
}

// specialisation for double
template<>
inline void RefCntStorage<double>::operator delete(void* stor) {
  freeRefCntStorageDouble(stor);
}

// specialisation for float
template<>
inline void* RefCntStorage<float>::operator new(size_t) {
  return allocRefCntStorageFloat();
}

// specialisation for float
template<>
inline void RefCntStorage<float>::operator delete(void* stor) {
  freeRefCntStorageFloat(stor);
}

template<typename T>
inline void RefCntStorage<T>::setExisting(const T* stor, T_Index len) {
  startAddr = const_cast<T*>( stor );
  size = len;
  refs = std::numeric_limits<T_Index>::max();
}

GLINEAR_NAMESPACE_END
