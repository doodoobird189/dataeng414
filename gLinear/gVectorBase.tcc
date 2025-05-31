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

// glinear headers
#include "gLinearMacros.h"

// standard headers
#include <iostream>  // cout, ostream, istream
#include <string>  // string


GLINEAR_NAMESPACE_BEGIN


template<typename T_Deriv>
const char gVectorBase<T_Deriv>::className[] = "gVectorBase<T_Vector>";



// dump info
template<typename T_Vector>
void gVectorBase<T_Vector>::info(int i) const {
  GLINEAR_STATIC_VIRTUAL(const T_Vector) specialInfo(i);
}

// return size
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Index
gVectorBase<T_Vector>::size(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialSize();
}

// resize
template<typename T_Vector>
bool gVectorBase<T_Vector>::resize(T_Index elem) {
  // return if resize unecessary
  if (size() == elem) {
    return false;
  } // if
  // now always check to catch misuse #ifdef CHK_IDX
  chkidxAlias("resize(T_Index rows, T_Index cols)", "resize");
  //#endif // CHK_IDX
  // update size
  GLINEAR_STATIC_VIRTUAL(T_Vector) specialResize(elem);
  return true;
}

// const transpose
template<typename T_Vector>
const typename gVectorBase<T_Vector>::T_Transpose
gVectorBase<T_Vector>::transpose() const {
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialTranspose();
}

// transpose
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Transpose
gVectorBase<T_Vector>::transpose() {
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialTranspose();
}

// vector assignment
template<typename T_Vector>
template<typename T_SrcVector>
T_Vector&
gVectorBase<T_Vector>::assignTo(const gVectorBase<T_SrcVector>& vec) {
  // catch self-assign
  if ( static_cast<const void*>(&vec) == static_cast<const void*>(this) ) {
    return reinterpret_cast<T_Vector&>(*this);
  } // if
  typedef typename gVectorBase<T_SrcVector>::T_ConstIterator T_SrcIterator;
  typedef typename T_Vector::T_Value tval;
  // attempt resize. (should catch any naughty behaviour).
  resize( vec.size() );
  T_SrcIterator src = vec.begin();
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ = static_cast<tval>( *src++ );
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// scalar assignment
template<typename T_Vector>
T_Vector&
gVectorBase<T_Vector>::assignToAll(T_ConstReference scalar) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ = scalar;
  } // while
  return reinterpret_cast<T_Vector&>(*this);
} // assignTo

// assign to set
template<typename T_Vector>
T_Vector&
gVectorBase<T_Vector>::set(T_ConstReference beginval,
         T_ConstReference increment) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  // make this an integer so we don't compound round-off error.
  T_Index inx=0;
  while (dst != dst_end) {
    *dst++ = beginval + (increment * static_cast<T_ConstReference>(inx++));
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// elem access
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Reference
gVectorBase<T_Vector>::operator[](T_Index elem) {
#ifdef CHK_IDX
  chkidxIndex("operator[](T_Index elem)", "Invalid element index", elem);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialElem(elem);
}

// const elem access
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_ConstReference
gVectorBase<T_Vector>::operator[](T_Index elem) const {
#ifdef CHK_IDX
  chkidxIndex("operator[](T_Index elem) const", "Invalid element index", elem);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialElem(elem);
}

// assign from C/C++ array
template<typename T_Vector>
void gVectorBase<T_Vector>::assignFromArray(T_ConstPointer addr,
              T_Index elem,
              T_Index stride) {
  if (elem == -1) {
    elem = size();
  } // if
  else {
    resize(elem);
  } // else
  T_Iterator dst = begin();
  for (T_Index indx = 0; indx < elem; indx++) {
    *dst++ = addr[indx*stride];
  } // for
}

// const begin iterator
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_ConstIterator
gVectorBase<T_Vector>::begin(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialBegin();
}

// begin iterator
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Iterator
gVectorBase<T_Vector>::begin(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialBegin();
}

// const past-the-end iterator
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_ConstIterator
gVectorBase<T_Vector>::end(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialEnd();
}

// past-the-end iterator
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Iterator
gVectorBase<T_Vector>::end(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialEnd();
}

// return const slice
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_ConstSlice
gVectorBase<T_Vector>::slice(const gIndexRange& elems) const {
#ifdef CHK_IDX
  chkidxSlice("gVectorBase<T_Vector>::slice(const gIndexRange& elems) const", "Invalid slice", elems);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialSlice(elems);
}

// return slice
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Slice
gVectorBase<T_Vector>::slice(const gIndexRange& elems) {
#ifdef CHK_IDX
  chkidxSlice("gVectorBase<T_Vector>::slice(const gIndexRange& elems)", "Invalid slice", elems);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialSlice(elems);
}

// return const subrange
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_ConstSubrange
gVectorBase<T_Vector>::subrange(T_Index fromElem,
        T_Index toElem) const {
#ifdef CHK_IDX
  chkidxSubrange("gVectorBase<T_Vector>::subrange(T_Index fromElem, T_Index toElem) const", "Invalid subrange",
     fromElem, toElem);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialSubrange(fromElem, toElem);
}

// return subrange
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Subrange
gVectorBase<T_Vector>::subrange(T_Index fromElem,
        T_Index toElem) {
#ifdef CHK_IDX
  chkidxSubrange("gVectorBase<T_Vector>::subrange(T_Index fromElem, T_Index toElem)", "Invalid subrange",
     fromElem, toElem);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Vector) specialSubrange(fromElem, toElem);
}

// return number of references
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Index
gVectorBase<T_Vector>::getRefs(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Vector) specialGetRefs();
}

// deepen
template<typename T_Vector>
inline void gVectorBase<T_Vector>::deepen(void) const {
  if (getRefs() == 1) {
    return;
  } // if
  GLINEAR_STATIC_VIRTUAL(const T_Vector) specialDeepen();
}

// alias
template<typename T_Vector>
inline void gVectorBase<T_Vector>::alias(T_Vector& vec) {
  GLINEAR_STATIC_VIRTUAL(T_Vector) specialAlias(vec);
}

// const alias
template<typename T_Vector>
inline void gVectorBase<T_Vector>::alias(const T_Vector& vec) const {
  GLINEAR_STATIC_VIRTUAL(const T_Vector) specialAlias(vec);
}

// minimum element
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::min(void) const {
#ifdef CHK_IDX
  chkidxZerosize("min(void) const", "Size should exceed 0");
#endif // CHK_IDX
  T_ConstIterator ptr = begin();
  const T_ConstIterator ptr_end = end();
  T_Value smallest = *ptr++;
  while (ptr != ptr_end) {
    if (*ptr < smallest) {
      smallest = *ptr;
    } // if
    ptr++;
  } // while
  return smallest;
}

// minimum element value in range
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::minVal(T_Index start,
            T_Index _end) const {
  if (_end == -1) {
    _end = size() - 1;  // default to end
  } // if
#ifdef CHK_IDX
  chkidxSubrange("minVal(T_Index start, T_Index end) const", "", start, _end);
#endif // CHK_IDX
  T_ConstIterator ptr = begin() + start;
  T_Value smallest = *ptr++;    // point to start + 1
  for (T_Index idx = start + 1; idx <= _end; idx++, ptr++) {
    if (*ptr < smallest) {
      smallest = *ptr;
    } // if
  } // for
  return smallest;
}

// minimum element position in range
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Index
gVectorBase<T_Vector>::minPos(T_Index start,
            T_Index _end) const {
  if (_end == -1) {
    _end = size() - 1;  // default to end
  } // if
#ifdef CHK_IDX
  chkidxSubrange("minPos(T_Index start, T_Index end) const", "", start, _end);
#endif // CHK_IDX
  T_Index smallestIdx = start;
  T_ConstIterator ptr = begin() + start;
  T_Value smallest = *ptr++;    // point to start + 1
  for (T_Index idx = start + 1; idx <= _end; idx++, ptr++) {
    if (*ptr < smallest) {
      smallestIdx = idx;
      smallest = *ptr;
    } // if
  } // for
  return smallestIdx;
}

// maximum element
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::max(void) const {
#ifdef CHK_IDX
  chkidxZerosize("max(void) const", "Size should exceed 0");
#endif // CHK_IDX
  T_ConstIterator ptr = begin();
  const T_ConstIterator ptr_end = end();
  T_Value largest = *ptr++;
  while (ptr != ptr_end) {
    if (*ptr > largest) {
      largest = *ptr;
    } // if
    ptr++;
  } // while
  return largest;
}

// maximum element value in range
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::maxVal(T_Index start,
            T_Index _end) const {
  if (_end == -1) {
    _end = size() - 1;  // default to end
  } // if
#ifdef CHK_IDX
  chkidxSubrange("maxVal(T_Index start, T_Index end) const", "", start, _end);
#endif // CHK_IDX
  T_ConstIterator ptr = begin() + start;
  T_Value largest = *ptr++;    // point to start + 1
  for (T_Index idx = start + 1; idx <= _end; idx++, ptr++) {
    if (*ptr > largest) {
      largest = *ptr;
    } // if
  } // for
  return largest;
}

// maximum element position in range
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Index
gVectorBase<T_Vector>::maxPos(T_Index start,
            T_Index _end) const {
  if (_end == -1) {
    _end = size() - 1;  // default to end
  } // if
#ifdef CHK_IDX
  chkidxSubrange("maxPos(T_Index start, T_Index end) const", "", start, _end);
#endif // CHK_IDX
  T_Index largestIdx = start;
  T_ConstIterator ptr = begin() + start;
  T_Value largest = *ptr++;    // point to start + 1
  for (T_Index idx = start + 1; idx <= _end; idx++, ptr++) {
    if (*ptr > largest) {
      largestIdx = idx;
      largest = *ptr;
    } // if
  } // for
  return largestIdx;
}

// maximum element
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::absMax(void) const {
#ifdef CHK_IDX
  chkidxZerosize("max(void) const", "Size should exceed 0");
#endif // CHK_IDX
  T_ConstIterator ptr = begin();
  const T_ConstIterator ptr_end = end();
  T_Value largest = *ptr++;
  if (largest < 0) {
    largest = -largest;
  } // if
  while (ptr != ptr_end) {
    if (*ptr > largest) {
      largest = *ptr;
    } // if
    if (-(*ptr) > largest) {
      largest = -(*ptr);
    } // if
    ptr++;
  } // while
  return largest;
}

// find elem
template<typename T_Vector>
typename gVectorBase<T_Vector>::T_Index
gVectorBase<T_Vector>::findElem(T_ConstReference wanted,
        T_Index startAt,
        int direction) const {
#ifdef CHK_IDX
  chkidxSubrange("findElem(const T_Value& wanted, T_Index startAt, int direction) const",
     "", startAt, size() - 1);
#endif // CHK_IDX
  // forward
  if (direction >= 0) {
    for (T_Index indx = startAt; indx < size(); indx++) {
      if (this->operator[](indx) == wanted) {
  return indx;
      } // if
    } // for
  } // if
  else { // reverse
    // OK, this is messy. But think it through while you start hacking away.
    // TIP: an unsigned can never be < 0.
    for (T_Index indx = startAt+1; indx-- > 0; ) {
      if (this->operator[](indx) == wanted) {
  return indx;
      } // if
    } // for
  } // else
  // not found
  return -1;
}

// vector computed assignment: addition
template<typename T_Vector>
template<typename T_SrcVector>
inline T_Vector&
gVectorBase<T_Vector>::operator+=(const gVectorBase<T_SrcVector>& vec) {
#ifdef CHK_IDX
  chkidxSize( "operator+=(const gVectorBase<T_SrcVector>& vec)",
        "", size(), vec.size() );
#endif // CHK_IDX
  typedef typename gVectorBase<T_SrcVector>::T_ConstIterator T_SrcIterator;
  T_SrcIterator src = vec.begin();
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *(dst++) += *(src++);
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// vector computed assignment: subtraction
template<typename T_Vector>
template<typename T_SrcVector>
inline T_Vector&
gVectorBase<T_Vector>::operator-=(const gVectorBase<T_SrcVector>& vec) {
#ifdef CHK_IDX
  chkidxSize( "operator-=(const gVectorBase<T_SrcVector>& vec)",
        "", size(), vec.size() );
#endif // CHK_IDX
  typedef typename gVectorBase<T_SrcVector>::T_ConstIterator T_SrcIterator;
  typedef typename gVectorBase<T_Vector>::T_Value T_DstValue;
  T_SrcIterator src = vec.begin();
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ -= static_cast<T_DstValue>(*src++);
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// scalar computed assignment: addition
template<typename T_Vector>
inline T_Vector&
gVectorBase<T_Vector>::operator+=(T_ConstReference scalar) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ += scalar;
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// scalar computed assignment: subtraction
template<typename T_Vector>
inline T_Vector&
gVectorBase<T_Vector>::operator-=(T_ConstReference scalar) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ -= scalar;
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// scalar computed assignment: multiplication
template<typename T_Vector>
inline T_Vector&
gVectorBase<T_Vector>::operator*=(T_ConstReference scalar) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ *= scalar;
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// scalar computed assignment: division
template<typename T_Vector>
inline T_Vector&
gVectorBase<T_Vector>::operator/=(T_ConstReference scalar) {
  T_Iterator dst = begin();
  const T_ConstIterator dst_end = end();
  while (dst != dst_end) {
    *dst++ /= scalar;
  } // while
  return reinterpret_cast<T_Vector&>(*this);
}

// sum
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::sum(void) const {
#ifdef CHK_IDX
  chkidxZerosize("sum(void) const", "Size should exceed 0");
#endif // CHK_IDX
  T_ConstIterator src = begin();
  const T_ConstIterator src_end = end();
  T_Value acc = *src++;
  while (src != src_end) {
    acc += *src++;
  } // while
  return acc;
}

// sum of squares
template<typename T_Vector>
inline typename gVectorBase<T_Vector>::T_Value
gVectorBase<T_Vector>::sum2(void) const {
#ifdef CHK_IDX
  chkidxZerosize("sum2(void) const", "Size should exceed 0");
#endif // CHK_IDX
  T_ConstIterator src = begin();
  const T_ConstIterator src_end = end();
  T_Value acc = (*src * *src);
  while (++src != src_end) {
    acc += (*src * *src);
  } // while
  return acc;
}

// dump to cout
template<typename T_Vector>
inline void gVectorBase<T_Vector>::print(int depth) const {
  if (depth < 0) {
    return;
  } // if
  std::cout << reinterpret_cast<const T_Vector&>(*this);
  /*
  int sz = size();
  for (int i =0; i < sz; i++) {
    std::cout << operator[](i) << ' ';
  }
  std::cout << std::endl;
  */
}


// binary output to stream
template<typename T_Vector>
std::ostream& gVectorBase<T_Vector>::write(std::ostream& os) const {
  return global_write( os, reinterpret_cast<const T_Vector&>(*this) );
}

// binary input from stream
template<typename T_Vector>
std::istream& gVectorBase<T_Vector>::read(std::istream& is) {
  return global_read( is, reinterpret_cast<T_Vector&>(*this) );
}

// index bound checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxIndex(std::string funcName,
          std::string infostr,
          T_Index elem) const {
  if ( (elem < 0) || ( elem >= size() ) ) {
    gErrorIndex(className, funcName, infostr, size()-1, elem);
  } // if
}

// zero size checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxZerosize(std::string funcName,
             std::string infostr) const {
  if (size() == 0) {
    gErrorSize( className, funcName, infostr, size(), size() );
  } // if
}

// size checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxSize(std::string funcName,
               std::string infostr,
               T_Index size1,
               T_Index size2) const {
  if (size1 != size2) {
    gErrorSize(className, funcName, infostr, size1, size2);
  } // if
}

// subrange index checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxSubrange(std::string funcName,
             std::string infostr,
             T_Index from,
             T_Index to) const {
  if ( (from > to) || ( to >= size() ) ) {
    gErrorSubrange(className, funcName, infostr, from, to, size() - 1);
  } // if
}

// slice index checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxSlice(std::string funcName,
          std::string infostr,
          const gIndexRange& elems) const {
  if ( ( elems.getStart() > elems.getEnd() ) || (elems.getEnd() > size() - 1) ) {
    gErrorSlice( className, funcName, infostr, size() - 1, elems.getStart(), elems.getEnd(), elems.getStride() );
  } // if
}
// alias checker
template<typename T_Vector>
void gVectorBase<T_Vector>::chkidxAlias(std::string funcName,
          std::string infostr) const {
  if (getRefs() > 1) {
    gErrorAliases(className, funcName, infostr);
  } // if
}


GLINEAR_NAMESPACE_END
