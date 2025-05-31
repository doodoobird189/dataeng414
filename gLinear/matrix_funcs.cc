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

char MATRIX_FUNCS_CC_VERSION[] = "$Id$";

// glinear headers
#include "matrix_funcs.h"
#include "exceptions.h"

// standard headers
#include <sstream>  // ostringstream
#include <iostream>  // endl
#ifdef ABORT_NOT_THROW
#include <cstdlib>  // abort
#endif // ABORT_NOT_THROW

GLINEAR_NAMESPACE_BEGIN

// isstolen.from("JA du Preez") == true;
/* inv:
 Gebruik LU dekomposisie
 */
gMatrix<float> inv(const gMatrix<float>& mat1,
                   float& invdeterm) {
   int nro = mat1.rows();
   int nco = mat1.cols();

   if (nro != nco) {
     std::ostringstream errMsg;
     errMsg << '\n' << "ERROR => " << "inv(const NumgMatrix<float>&) : matrix not square "
            << nro << 'x' << nco << std::endl;
     errMsg << __FILE__ << ":" << __LINE__ << std::endl;
#ifdef ABORT_NOT_THROW
     std::cerr << errMsg.str() << std::flush;
     abort();
#else
     throw GLinearError( errMsg.str() );
#endif // ABORT_NOT_THROW
   } // if

   gMatrix<float> antw(nro, nco), mat2(nro,nco);
   mat2 = mat1;
   gRowVector<int> indx(nro);
   float d(1.0);
   gRowVector<float> col(nro);
   ludcmp(mat2, nro, indx, d);

   invdeterm = mat2(0, 0);
   for (int i = 1; i < nro; i++) {
     invdeterm *= mat2(i, i);
   } // for i
   invdeterm *= d;
   invdeterm = 1.0f / invdeterm;

   for (int j = 0; j < nro; j++) {
     for (int i = 0; i < nro; i++) {
       col[i] = 0;
     } // for i
     col[j] = 1;
     lubksb(mat2, nro, indx, col);
     for (int i = 0; i < nro; i++) {
       antw(i, j) = col[i];
     } // for i
   } // for j
   return antw;
} // inv


#define TINY 1.0e-20f;
// isstolen.from("JA du Preez") == true;
void ludcmp(gMatrix<float>& a,
      int n,
      gRowVector<int>& indx,
      float& d) {

  gRowVector<float> vv(n);
  d = 1.0;
  for (int i = 0; i < n; i++) {
    float big(0.0);
    for (int j = 0; j < n; j++) {
      float temp = fabsf( a(i, j) );
      if (temp > big) {
        big = temp;
      } // if
    } // for
    if (big == 0.0) {
      std::ostringstream errMsg;
      errMsg << '\n' << "ERROR => " << "Singular matrix in routine MATRIX :LUDCMP" << std::endl;
      errMsg << __FILE__ << ":" << __LINE__ << std::endl;
#ifdef ABORT_NOT_THROW
      std::cerr << errMsg.str() << std::flush;
      abort();
#else
      throw GLinearError( errMsg.str() );
#endif // ABORT_NOT_THROW
    } // if
    vv[i] = 1.0f / big;
  } // for
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      float sum = a(i, j);
      for (int k = 0; k < i; k++) {
        sum -= a(i, k) * a(k, j);
      } // for k
      a(i, j) = sum;
    } // for i
    float big(0.0);
    int imax(-1);
    for (int i = j; i < n; i++) {
      float sum = a(i, j);
      for (int k = 0; k < j; k++) {
        sum -= a(i, k) * a(k, j);
      } // for k
      a(i, j) = sum;
      float dum = vv[i] * fabsf(sum);
      if (dum >= big) {
        big = dum;
        imax = i;
      } // if
    }  // for i
    if (j != imax) {
      for (int k = 0; k < n; k++) {
        float dum = a(imax, k);
        a(imax, k) = a(j, k);
        a(j, k) = dum;
      } // for
      d = -(d);
      vv[imax] = vv[j];
    } // if
    indx[j] = imax;
    if (a(j, j) == 0.0) {
      a(j, j) = TINY;
    } // if
    if (j != n) {
      float dum = 1.0f / a(j, j);
      for (int i = j + 1; i < n; i++) {
        a(i, j) *= dum;
      } // for
    } // if
  } // for j
} // ludcmp

#undef TINY

void lubksb(gMatrix<float>& a,
      int n,
      gRowVector<int>& indx,
      gRowVector<float>& b) {
  int ii(-1);

  for (int i = 0; i < n; i++) {
    int ip = indx[i];
    float sum = b[ip];
    b[ip] = b[i];
    if (ii + 1)
      for (int j = ii; j <= i - 1; j++) {
        sum -= a(i, j) * b[j];
      } // for j
    else if (sum) {
      ii = i;
    } // else if
    b[i] = sum;
  } // for i
  for (int i = n - 1; i >= 0; i--) {
    float sum = b[i];
    for (int j = i + 1; j < n; j++) {
      sum -= a(i, j) * b[j];
    } // for
    b[i] = sum / a(i, i);
  } // for i
} // lubksb

GLINEAR_NAMESPACE_END
