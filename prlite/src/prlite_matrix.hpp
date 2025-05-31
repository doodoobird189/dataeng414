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


#ifndef PRLITE_MATRIX_HPP
#define PRLITE_MATRIX_HPP

// patrec headers
#include "prlite_gendefs.hpp"  // REAL
#include "prlite_genmat.hpp"  // Matrix, DiagMatrix
#include "prlite_genvec.hpp"  // Vector

// glinear headers
#include "gLinear/gMatrixOpers.h"
// #include "../gLinear/matrix_funcs.h"

namespace prlite{

template<typename T>
T calcDet(const gLinear::gRowMatrix<T>& mat1);

extern gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<float>& v);
extern gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<double>& v);
extern gLinear::gRowMatrix<REAL> outerProd(const gLinear::gRowVector<float>& cvec,
                              const gLinear::gRowVector<float>& rvec);
extern void outerProd(gLinear::gRowMatrix<REAL>& outer,
                      const gLinear::gRowVector<float>& cvec,
                      const gLinear::gRowVector<float>& rvec);

extern prlite::DiagMatrix<REAL> inv(const prlite::DiagMatrix<REAL>& mat1,
                            REAL& invdet);

/*
extern REAL quadraticForm(const Vector<REAL>& x,
                          const Matrix<REAL>& m);
*/
/**
 * Calculates
 * \f$\sum_i (x[i])^2m[i,i] + 2\sum_i\sum_{j>i}x[i]x[j]m[i,j]\f$
 */
/*
extern REAL quadraticForm(const Vector<REAL>& x,
                          const prlite::DiagMatrix<REAL>& m);

extern Vector<float> operator*(const Vector<float>& v,
                               const prlite::Transpose<float>& m); // row*mat
extern void MPY(Vector<float>& prj,
                const Vector<float>& v,
                const prlite::Transpose<float>& m); // row*mat
*/

template<typename T>
gLinear::gRowVector<T> asVector(const gLinear::gRowMatrix<T>& mat);

} // namespace prlite

#include "prlite_matrix.tcc"


#endif // PRLITE_MATRIX_HPP
