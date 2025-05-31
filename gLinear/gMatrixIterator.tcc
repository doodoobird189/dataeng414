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

template<typename T_Matrix>
inline typename gMatrixIteratorSupport<T_Matrix, ROW>::T_Row
gMatrixIteratorSupport<T_Matrix, ROW>::deref(const T_Iter& iter) {
  return iter.matrix->row(iter.indx);
}


template<typename T_Matrix>
inline typename gMatrixIteratorSupport<T_Matrix, COL>::T_Col
gMatrixIteratorSupport<T_Matrix, COL>::deref(const T_Iter& iter) {
  return iter.matrix->col(iter.indx);
}

// default ctor
template<typename T_Matrix, typename T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::gMatrixIterator(void) :
  indx(0),
  matrix(NULL) {
}

// init ctor
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>::gMatrixIterator(
  T_Index new_indx, const T_Matrix* new_matrix) :
  indx(new_indx),
  matrix( const_cast<T_Matrix*>(new_matrix) ) {
}

// copy ctor
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>::gMatrixIterator(
  const gMatrixIterator<T_Matrix, T_VectorTag>& iter) :
  indx(iter.indx),
  matrix(iter.matrix) {
}
#if 0
// const copy assignment
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator=(const gMatrixIterator<T_Matrix, T_VectorTag>& iter) const {
  indx = iter.indx;
  matrix = iter.matrix;
  return *this;
}
# endif
// copy assignment
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator=(gMatrixIterator<T_Matrix, T_VectorTag>& iter) {
  indx = iter.indx;
  matrix = iter.matrix;
  return *this;
}

// const dereference
template<typename T_Matrix, typename T_VectorTag>
const typename gMatrixIterator<T_Matrix, T_VectorTag>::T_Vector
gMatrixIterator<T_Matrix, T_VectorTag>::operator*(void) const {
  return T_Support::deref(*this);
}

// dereference
template<typename T_Matrix, typename T_VectorTag>
typename gMatrixIterator<T_Matrix, T_VectorTag>::T_Vector
gMatrixIterator<T_Matrix, T_VectorTag>::operator*(void) {
  return gMatrixIteratorSupport<T_Matrix, T_VectorTag>::deref(*this);
}

// equality operator
template<typename T_Matrix, typename T_VectorTag>
inline bool
gMatrixIterator<T_Matrix, T_VectorTag>::operator==(const T_This& iter) {
  return (indx == iter.indx);
}

// inequality operator
template<typename T_Matrix, typename T_VectorTag>
inline bool
gMatrixIterator<T_Matrix, T_VectorTag>::operator!=(const T_This& iter) {
  return !(indx == iter.indx);
}

// postfix increment support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::postfixInc(void) {
  T_This tmp(*this);
  indx++;
  return tmp;
}

// prefix increment support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::prefixInc(void) {
  indx++;
  return *this;
}

// postfix decrement support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::postfixDec(void) {
  T_This tmp(*this);
  indx--;
  return tmp;
}

// prefix decrement support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::prefixDec(void) {
  indx--;
  return *this;
}
#if 0
// const postfix increment
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator++(int) const {
  return postfixInc();
}
#endif
// postfix increment
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator++(int) {
  return postfixInc();
}
#if 0
// const prefix increment
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator++(void) const {
  return prefixInc();
}
#endif
// prefix increment
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator++(void) {
  return prefixInc();
}
#if 0
// const postfix decrement
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator--(int) const {
  return postfixDec();
}
#endif
// postfix decrement
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator--(int) {
  return postfixDec();
}
#if 0
// const prefix decrement
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator--(void) const {
  return prefixDec();
}
#endif
// prefix decrement
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>&
gMatrixIterator<T_Matrix, T_VectorTag>::operator--(void) {
  return prefixDec();
}

// addition support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::add(T_Index i) const {
  return T_This(indx + i, matrix);
}

// subtraction support
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::subtract(T_Index i) const {
  return T_This(indx - i, matrix);
}

// const addition
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator+(T_Index i) const {
  return add(i);
}

// addition
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator+(T_Index i) {
  return add(i);
}

// const subtraction
template<typename T_Matrix, typename T_VectorTag>
inline const gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator-(T_Index i) const {
  return subtract(i);
}

// subtraction
template<typename T_Matrix, typename T_VectorTag>
inline gMatrixIterator<T_Matrix, T_VectorTag>
gMatrixIterator<T_Matrix, T_VectorTag>::operator-(T_Index i) {
  return subtract(i);
}


GLINEAR_NAMESPACE_END
