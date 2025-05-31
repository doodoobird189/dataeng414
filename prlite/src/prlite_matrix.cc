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

char MATRIX_CC_VERSION[] = "$Id$";

/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_matrix.hpp"

using namespace std;

namespace prlite{

/******************************************************************************

IMPLEMENTATION NOTES:
The following routines has been optimized for speed by eliminating indexing
operators as far as practical.

*******************************************************************************/

gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<float>& v) {
  int sz = v.size();
  gLinear::gRowMatrix<REAL> outer(sz, sz);
  int i;
  for (i = 0; i < sz; i++) {
    for (int j = i; j < sz; j++) {
      outer[i][j] = v[i] * v[j];
    } // for j
  } // for i
  for (i = 1; i < sz; i++) {
    for (int j = 0; j < i; j++) {
      outer[i][j] = outer[j][i];
    } // for j
  } // for i
  return outer;
} // outerProd


gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<double>& v) {
  int sz = v.size();
  gLinear::gRowMatrix<REAL> outer(sz, sz);
  int i;
  for (i = 0; i < sz; i++) {
    for (int j = i; j < sz; j++) {
      outer[i][j] = v[i] * v[j];
    } // for j
  } // for i
  for (i = 1; i < sz; i++) {
    for (int j = 0; j < i; j++) {
      outer[i][j] = outer[j][i];
    } // for j
  } // for i
  return outer;
} // outerProd


gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<float>& cvec,
                       const gLinear::gRowVector<float>& rvec) {
  gLinear::gRowMatrix<REAL> outer( cvec.size(), rvec.size() );
  outerProd(outer, cvec, rvec);
  return outer;
} // outerProd


void outerProd(gLinear::gRowMatrix<REAL>& outer,
               const gLinear::gRowVector<float>& cvec,
               const gLinear::gRowVector<float>& rvec) {
  int nro = cvec.size();
  int nco = rvec.size();
  if (outer.rows() != nro || outer.cols() != nco) {
    PRLITE_FAIL("outerProd(Matrix<REAL>&, Vector<float>&, Vector<float>&) : arguments incompatible");
  } // if
  int oRowIdx = 0;
  // const float* cPtr = cvec.arrayAddress();
  gLinear::gRowVector<float>::T_ConstIterator cItr = cvec.begin();
  while (nro--) {
    gLinear::gRowVector<float>::T_ConstIterator rItr = rvec.begin();
    // const float* rPtr = rvec.arrayAddress();
    gLinear::gRowVector<REAL>::T_Iterator oRowItr = outer[oRowIdx].begin();
    // float* oRowPtr = outer[oRowIdx].arrayAddress();
    int j = nco;
    while (j--) {
      *oRowItr++ = *cItr * *rItr++;
    } // while j
    cItr++;
    oRowIdx++;
  } // while
} // outerProd


prlite::DiagMatrix<REAL> inv(const prlite::DiagMatrix<REAL>& mat1,
                     REAL& invdeterm) { // determinant at same time
  int nro = mat1.size();
  prlite::DiagMatrix<REAL> antw(nro);

  antw[0] = 1.0 / mat1[0];
  invdeterm = antw[0];
  for (int i = 1; i < nro; i++) {
    antw[i] = 1.0 / mat1[i];
    invdeterm *= antw[i];
  } // for i

  return antw;
} // inv

} // namespace prlite
