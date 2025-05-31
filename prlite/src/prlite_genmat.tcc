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

/*
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : pre-2000
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_genvec.hpp"  // Vector
#include "prlite_error.hpp"  // FAIL

// standard headers
#include <iostream>  // istream, ostream, endl, cout
#include <algorithm>  // min
#include <string>  // string

namespace prlite{

template<typename T>
inline int length(const RowMatrix<T>& mat) {
  return mat.rows();
} // length

template<typename T>
inline std::istream& operator>>(std::istream& file,
        DiagMatrix<T>& dm) {
  PRLITE_CHECK(file >> dm.mainDiag, IOError, "Error reading Vector.");
  return file;
} // operator>>

template<typename T>
inline std::ostream& operator<<(std::ostream& file,
        const DiagMatrix<T>& dm) {
  return file << dm.mainDiag << std::endl;
} // operator<<

template<typename T>
inline DiagMatrix<T>::DiagMatrix() {
  offDiag = 0;
  /*
    odPtr = new T;
    *odPtr = 0;
    */
} // default constructor

template<typename T>
inline DiagMatrix<T>::DiagMatrix(const DiagMatrix<T>& m)
    : mainDiag(m.mainDiag) {
  offDiag = 0;
  /*
    odPtr = new T;
    *odPtr = 0;
    */
} // copy constructor

template<typename T>
inline DiagMatrix<T>::~DiagMatrix() {
  /*
    delete odPtr;
    odPtr = 0;
  */
} // destructor

template<typename T>
DiagMatrix<T>& DiagMatrix<T>::operator=(const DiagMatrix<T>& m) {
  mainDiag = m.mainDiag;
  return *this;
} // operator=

template<typename T>
DiagMatrix<T>& DiagMatrix<T>::operator=(const T& m) {
  mainDiag = m;
  return *this;
} // operator=

template<typename T>
inline DiagMatrix<T>::DiagMatrix(int dim)
    : mainDiag(dim) {
  offDiag = 0;
  /*
    odPtr = new T;
    *odPtr = 0;
    */
} // constructor

template<typename T>
DiagMatrix<T>::DiagMatrix(const RowMatrix<T>& m)
    : mainDiag( std::min( m.rows(), m.cols() ) ) {
  if ( m.rows() != m.cols() ) {
    std::cout << "WARNING => DiagMatrix(const Matrix&) : non-square matrix (" << m.rows() << 'x' << m.cols() << std::endl;
  } // if
  for (int i = 0; i < mainDiag.size(); i++) {
    mainDiag[i] = m[i][i];
  } // for i
  offDiag = 0;
  /*
    odPtr = new T;
    *odPtr = 0;
    */
} // constructor

template<typename T>
inline const std::string& DiagMatrix<T>::isA() const {
  static const std::string CLASSNAME("DiagMatrix<T>");
  return CLASSNAME;
} // isA

template<typename T>
inline void DiagMatrix<T>::resize(int sz) {
  mainDiag.resize(sz);
} // resize

template<typename T>
inline void DiagMatrix<T>::deepen() {
  mainDiag.deepen();
} // deepen

template<typename T>
inline void DiagMatrix<T>::setAll(const T& val) {
  mainDiag = (val);
} // setAll

template<typename T>
inline T& DiagMatrix<T>::operator[](int i) {
  return mainDiag[i];
} // operator[]

template<typename T>
inline T DiagMatrix<T>::operator[](int i) const {
  return mainDiag[i];
} // operator[]

template<typename T>
inline T DiagMatrix<T>::operator()(int i) const {
  return mainDiag[i];
} // operator()

template<typename T>
inline T& DiagMatrix<T>::operator()(int i) {
  return mainDiag[i];
} // operator()

template<typename T>
inline T DiagMatrix<T>::elem(int i) const {
  return mainDiag[i];
} // elem

template<typename T>
inline T& DiagMatrix<T>::elem(int i) {
  return mainDiag[i];
} // elem

template<typename T>
DiagMatrix<T>::operator RowMatrix<T>() const {
  int sz( size() );
  RowMatrix<T> tmp(sz, sz);
  tmp.assignToAll(0); // .assignToAll(0); // setAll(0.0F) ??
  for (int i = 0; i < sz; i++) {
    tmp[i][i] = mainDiag[i];
  } // for i
  return tmp;
} // operator Matrix<T>

template<typename T>
int DiagMatrix<T>::size() const {
  return mainDiag.size();
} // size

template<typename T>
const T* DiagMatrix<T>::arrayAddress() const {
  return mainDiag.arrayAddress();
} // arrayAddress

template<typename T>
T* DiagMatrix<T>::arrayAddress() {
  return mainDiag.arrayAddress();
} // arrayAddress

template<typename T>
T DiagMatrix<T>::operator()(int r,
                            int c) const {
  if (r == c) {
    return mainDiag[r];
  } // if
  else {
    offDiag = 0;
    return offDiag;
    /*
     *odPtr = 0;
     return *odPtr;
    */
  } // else
} // operator()

template<typename T>
T& DiagMatrix<T>::operator()(int r,
                             int c) {
  if (r == c) {
    return mainDiag[r];
  } // if
  else {
    offDiag = 0;
    return offDiag;
    /*
     *odPtr = 0;
     return *odPtr;
    */
  } // else
} // operator()

template<typename T>
T DiagMatrix<T>::elem(int r,
                      int c) const {
  if (r == c) {
    return mainDiag[r];
  } // if
  else {
    offDiag = 0;
    return offDiag;
    /*
     *odPtr = 0;
     return *odPtr;
    */
  } // else
} // elem

template<typename T>
T& DiagMatrix<T>::elem(int r,
                       int c) {
  if (r == c) {
    return mainDiag[r];
  } // if
  else {
    offDiag = 0;
    return offDiag;
    /*
     *odPtr = 0;
     return *odPtr;
    */
  } // else
} // elem

template<typename T>
void DiagMatrix<T>::print(int depth) const {
  if (depth < 0) {
    return;
  } // if
  mainDiag.print(depth);
} // print

template<typename T>
std::ostream& DiagMatrix<T>::write(std::ostream& file) const {
  mainDiag.write(file);
  return file;
} // write

template<typename T>
std::istream& DiagMatrix<T>::read(std::istream& file) {
  mainDiag.read(file);
  return file;
} // read

template<typename T>
typename RowVector<T>::T_ConstIterator DiagMatrix<T>::begin() const {
  return mainDiag.begin();
} // begin

template<typename T>
typename RowVector<T>::T_ConstIterator DiagMatrix<T>::end() const {
  return mainDiag.end();
} // end

template<typename T>
typename RowVector<T>::T_Iterator DiagMatrix<T>::begin() {
  return mainDiag.begin();
} // begin

template<typename T>
typename RowVector<T>::T_Iterator DiagMatrix<T>::end() {
  return mainDiag.end();
} // end

template<typename T>
inline Transpose<T>::Transpose(const RowMatrix<T>& mat)
    : orig(mat) { // for efficiency - no determinant
} // constructor

template<typename T>
inline int Transpose<T>::rows() const {
  return orig.cols();
} // rows

template<typename T>
inline int Transpose<T>::cols() const {
  return orig.rows();
} // cols

template<typename T>
inline T Transpose<T>::elem(int i,
                            int j) const {
  return orig.elem(j, i);
} // elem

template<typename T>
inline T& Transpose<T>::elem(int i,
                             int j) {
  return orig.elem(j, i);
} // elem

template<typename T>
T Transpose<T>::operator()(int i,
                           int j) const {
  if (j < 0 || orig.rows() <= j) {
    PRLITE_FAIL(isA() << "::operator()(int, ->int<-)" << " index out of bounds " << orig.rows() << j);
  } // if
  if (i < 0 || orig.cols() <= i) {
    PRLITE_FAIL(isA() << "::operator()(->int<-, int)" << " index out of bounds " << orig.cols() << i);
  } // if
  return orig.operator[](i)[j];
} // operator()

template<typename T>
T& Transpose<T>::operator()(int i,
                            int j) {
  if (j < 0 || orig.rows() <= j) {
    PRLITE_FAIL(isA() << "::operator()(int, ->int<-)" << " index out of bounds " << orig.rows() << j);
  } // if
  if (i < 0 || orig.cols() <= i) {
    PRLITE_FAIL(isA() << "::operator()(->int<-, int)" << " index out of bounds " << orig.cols() << i);
  } // if
  return orig.operator[](i)[j];
} // operator()

template<typename T>
Transpose<T>::operator const RowMatrix<T>() const {
  int co( orig.rows() );
  int ro( orig.cols() );
  RowMatrix<T> mat(ro, co);
  for (int i = 0; i < ro; i++) {
    for (int j = 0; j < co; j++) {
      mat(i, j) = orig(j, i);
    } // for j
  } // for i
  return mat;
} // operator const Matrix<T>

/*
template<typename T>
Vector<T>* Transpose<T>::arrayAddress() const {
  return (Vector<T>*) orig.arrayAddress();
} // arrayAddress
*/

template<typename T>
const std::string& Transpose<T>::isA() const {
  static const std::string CLASSNAME("Transpose<T>");
  return CLASSNAME;
} // isA

} // namespace prlite

/*
template<typename T>
Transpose<T>::operator const Matrix<T>() const {
  int co( orig.rows() );
  int ro( orig.cols() );
  Matrix<T> mat(ro, co);
  Vector<T>* basePtr( mat.arrayAddress() );
  for (int i = 0; i < ro; i++) {
    T* rowPtr( basePtr->arrayAddress() );
    basePtr++;
    for (int j = 0; j < co; j++) {
      *rowPtr++ = orig.operator[](j)[i];
    } // for j
  } // for i
  return mat;
} // operator const Matrix<T>
*/
