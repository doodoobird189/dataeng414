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
#include "prlite_genmat.hpp"  // Matrix, det
#include "prlite_genvec.hpp"  // Vector

namespace prlite{

template<typename T>
T calcDet(const gLinear::gRowMatrix<T>& mat1) {
  int fail;
  T theDet( det(mat1, fail) );
  PRLITE_ASSERT(!fail, "Matrix determinant calculation failed.");
  return theDet;
} // calcDet

template<typename T>
gLinear::gRowVector<T> asVector(const gLinear::gRowMatrix<T>& mat) {

  gLinear::gRowVector<T> tmp( mat.rows() * mat.cols() );
  int indx = 0;

  for (int row = 0; row < mat.rows(); row++) {
    for (int col = 0; col < mat.cols(); col++) {
      tmp[indx++] = mat(row, col);
    } // for (cols)
  } // for (rows)

  return tmp;
} // asVector

} // namespace prlite
