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


GLINEAR_NAMESPACE_BEGIN

// default ctor
template<typename T>
gSliceVectorIterator<T>::gSliceVectorIterator(void) :
  std::iterator<std::bidirectional_iterator_tag, T, T_Index>(),
  addr(NULL),
  stride(0) {
}

// ctor and init.
template<typename T>
gSliceVectorIterator<T>::gSliceVectorIterator(T_ConstPointer new_addr,
                                              T_Index new_stride) :
  std::iterator<std::bidirectional_iterator_tag, T, T_Index>(),
  addr( const_cast<T_Pointer>(new_addr) ),
  stride(new_stride) {
}

// copy ctor.
template<typename T>
gSliceVectorIterator<T>::gSliceVectorIterator(const gSliceVectorIterator<T>& iter) :
  std::iterator<std::bidirectional_iterator_tag, T, T_Index>(iter),
  addr(iter.addr),
  stride(iter.stride) {
}
// assigment operator
template<typename T>
gSliceVectorIterator<T>
gSliceVectorIterator<T>::operator=(const gSliceVectorIterator<T>& iter) {
  addr = iter.addr;
  stride = iter.stride;
  return *this;
}

// const dereference
template<typename T>
typename gSliceVectorIterator<T>::T_ConstReference
gSliceVectorIterator<T>::operator*(void) const {
  return *addr;
}

// dereference
template<typename T>
typename gSliceVectorIterator<T>::T_Reference
gSliceVectorIterator<T>::operator*(void) {
  return *addr;
}

// equality operator
template<typename T>
bool gSliceVectorIterator<T>::operator==(const gSliceVectorIterator<T>& iter) const {
  return (addr == iter.addr);
}

// inequality operator
template<typename T>
bool gSliceVectorIterator<T>::operator!=(const gSliceVectorIterator<T>& iter) const {
  return (addr != iter.addr);
}

// postfix increment
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator++(int) {
  gSliceVectorIterator<T> tmp(*this);
  addr += stride;
  return tmp;
}

// prefix increment
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator++(void) {
  addr += stride;
  return *this;
}

// postfix decrement
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator--(int) {
  gSliceVectorIterator<T> tmp(*this);
  addr -= stride;
  return tmp;
}

// prefix decrement
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator--(void) {
  addr -= stride;
  return *this;
}

// const addition
template<typename T>
const gSliceVectorIterator<T> gSliceVectorIterator<T>::operator+(T_Index delta) const {
  return gSliceVectorIterator<T>(addr + delta * stride, stride);
}

// addition
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator+(T_Index delta) {
  return gSliceVectorIterator<T>(addr + delta * stride, stride);
}

// const subtraction
template<typename T>
const gSliceVectorIterator<T> gSliceVectorIterator<T>::operator-(T_Index delta) const {
  return gSliceVectorIterator<T>(addr - delta * stride, stride);
}

// subtraction
template<typename T>
gSliceVectorIterator<T> gSliceVectorIterator<T>::operator-(T_Index delta) {
  return gSliceVectorIterator<T>(addr - delta * stride, stride);
}

GLINEAR_NAMESPACE_END
