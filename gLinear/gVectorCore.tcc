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
#include <iostream>  // cout, endl
#include <cmath>  // ceil
#include <new>  // placement new


GLINEAR_NAMESPACE_BEGIN

template<typename T_Vector>
const char gDenseVector<T_Vector>::className[] = "gDenseVector<T_Vector>";


template<typename T_Vector>
const char gSliceVector<T_Vector>::className[] = "gSliceVector<T_Vector>";


template<typename TF_Vector>
inline gDenseVector<TF_Vector>::gDenseVector(gLinearNoInit)
    : T_Base( gLinearNoInit() ) {
}

template<typename TF_Vector>
inline gDenseVector<TF_Vector>::gDenseVector(T_Index elem)
    : T_Base(elem) {
}

template<typename TF_Vector>
inline gDenseVector<TF_Vector>::gDenseVector(const gDenseVector<T_Vector>& vec) :
  T_Base(vec) {
}

template<typename TF_Vector>
inline gDenseVector<TF_Vector>::gDenseVector(T_Index new_nElem,
               const T_StorageElem* new_stor,
               T_ConstPointer new_startAddr) :
  T_Base(new_nElem, new_stor, new_startAddr) {
}

template<typename TF_Vector>
inline gDenseVector<TF_Vector>::gDenseVector(T_ConstPointer addr,
               T_Index elem) :
  T_Base(addr, elem, vector_no_copy_tag()) {
}

template<typename TF_Vector>
inline void gDenseVector<TF_Vector>::configure(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr) {
  this->~gDenseVector();
  new(this) gDenseVector(new_nElem, new_stor, new_startAddr);
}

template<typename TF_Vector>
inline void gDenseVector<TF_Vector>::specialInfo(int i) const {
  std::cout << "gDenseVector<T_Vector>" << std::endl;
  gVectorMemBase<T_Vector>::specialInfo(i);
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Reference
gDenseVector<TF_Vector>::specialElem(T_Index elem) const {
  return const_cast<T_Reference>( *(this->getStartAddr() + elem) );
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Iterator
gDenseVector<TF_Vector>::specialBegin(void) const {
  return const_cast<T_Iterator>( this->getStartAddr() );
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Iterator
gDenseVector<TF_Vector>::specialEnd(void) const {
  return const_cast<T_Iterator>( this->getStartAddr() + this->size() );
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Slice
gDenseVector<TF_Vector>::specialSlice(const gIndexRange& elems) const {
  T_Index new_nElem = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + elems.getEnd() - elems.getStart() ) /
                   static_cast<double>( elems.getStride() ) ) );
  T_ConstPointer new_startAddr = this->getStartAddr() + elems.getStart();
  return T_Slice( new_nElem, this->getStor(), new_startAddr, elems.getStride() );
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Subrange
gDenseVector<TF_Vector>::specialSubrange(T_Index fromElem,
           T_Index toElem) const {
  return T_Subrange(1 + toElem - fromElem, this->getStor(), this->getStartAddr() + fromElem);
}

template<typename TF_Vector>
inline typename gDenseVector<TF_Vector>::T_Transpose
gDenseVector<TF_Vector>::specialTranspose(void) const {
  return T_Transpose( this->getNElem(), this->getStor(), this->getStartAddr() );
}

template<typename TF_Vector>
inline gSliceVector<TF_Vector>::gSliceVector(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Vector>
inline gSliceVector<TF_Vector>::gSliceVector(T_Index elem) :
  T_Base(elem) {
  std::cout << "gSliceVector::gSliceVector() " << elem << std::endl;
}

template<typename TF_Vector>
inline gSliceVector<TF_Vector>::gSliceVector(const gSliceVector<T_Vector>& vec) :
  T_Base(vec),
  stride(vec.stride) {
}

template<typename TF_Vector>
inline gSliceVector<TF_Vector>::gSliceVector(T_Index new_nElem,
               const T_StorageElem* new_stor,
               T_ConstPointer new_startAddr,
               T_Index new_stride) :
  T_Base(new_nElem, new_stor, new_startAddr),
  stride(new_stride) {
}

template<typename TF_Vector>
inline void gSliceVector<TF_Vector>::configure(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr,
                 T_Index new_stride) {
  this->~gSliceVector();
  new(this) gSliceVector(new_nElem, new_stor, new_startAddr, new_stride);
}

template<typename TF_Vector>
inline void gSliceVector<TF_Vector>::specialInfo(int i) const {
  std::cout << "gSliceVector<T, gVector_SLICE>" << std::endl;
  std::cout << "  stride: " << stride << std::endl;
  gVectorMemBase<T_Vector>::specialInfo(i);
}

template<typename TF_Vector>
inline bool gSliceVector<TF_Vector>::specialResize(T_Index) {
  gErrorSliceResize(className, "resize(T_Index elem)");
  return true;
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Reference
gSliceVector<TF_Vector>::specialElem(T_Index elem) const {
  return const_cast<T_Reference>( this->getStartAddr()[elem * stride] );
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Iterator
gSliceVector<TF_Vector>::specialBegin(void) const {
  return gSliceVectorIterator<T_Value>( this->getStartAddr(), getStride() );
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Iterator
gSliceVector<TF_Vector>::specialEnd(void) const {
  return gSliceVectorIterator<T_Value>( this->getStartAddr(), getStride() ) + this->size();
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Slice
gSliceVector<TF_Vector>::specialSlice(const gIndexRange& elems) const {
  T_Index new_nElem = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + elems.getEnd() - elems.getStart() ) /
                   static_cast<double>( elems.getStride() ) ) );
  T_ConstPointer new_startAddr = this->getStartAddr() + ( stride * elems.getStart() );
  T_Index new_stride = stride * elems.getStride();
  return T_Slice(new_nElem, this->getStor(), new_startAddr, new_stride);
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Subrange
gSliceVector<TF_Vector>::specialSubrange(T_Index fromElem,
           T_Index toElem) const {
  return specialSlice( gIndexRange(fromElem, toElem, 1) );
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Transpose
gSliceVector<TF_Vector>::specialTranspose(void) const {
  return T_Transpose( this->getNElem(), this->getStor(), this->getStartAddr(), getStride() );
}

template<typename TF_Vector>
inline typename gSliceVector<TF_Vector>::T_Index
gSliceVector<TF_Vector>::getStride(void) const {
  return stride;
}

GLINEAR_NAMESPACE_END
