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

#ifndef GLAPACK_HEADERS_H
#define GLAPACK_HEADERS_H

// glinear headers
#include "gFortranTraits.h"

namespace fortran {

   extern "C" {
      // Eigens
      void cgeevx_(char *balanc,
                  char *jobvl,
                  char *jobvr,
                  char *sense,
                  integer *n,
                  complex *a,
                  integer *lda,
                  complex *w,
                  complex *vl,
                  integer *ldvl,
                  complex *vr,
                  integer *ldvr,
                  integer *ilo,
                  integer *ihi,
                  real *scale,
                  real *abnrm,
                  real *rconde,
                  real *rcondv,
                  complex *work,
                  integer *lwork,
                  real *rwork,
                  integer *info);

      void zgeevx_(char *balanc,
                  char *jobvl,
                  char *jobvr,
                  char *sense,
                  integer *n,
                  doublecomplex *a,
                  integer *lda,
                  doublecomplex *w,
                  doublecomplex *vl,
                  integer *ldvl,
                  doublecomplex *vr,
                  integer *ldvr,
                  integer *ilo,
                  integer *ihi,
                  doublereal *scale,
                  doublereal *abnrm,
                  doublereal *rconde,
                  doublereal *rcondv,
                  doublecomplex *work,
                  integer *lwork,
                  doublereal *rwork,
                  integer *info);

      // SVD
      void sgesvd_(char *jobu,
                  char *jobvt,
                  integer *m,
                  integer *n,
                  real *a,
                  integer *lda,
                  real *s,
                  real *u,
                  integer *ldu,
                  real *vt,
                  integer *ldvt,
                  real *work,
                  integer *lwork,
                  integer *info);

      void dgesvd_(char *jobu,
                  char *jobvt,
                  integer *m,
                  integer *n,
                  doublereal *a,
                  integer *lda,
                  doublereal *s,
                  doublereal *u,
                  integer *ldu,
                  doublereal *vt,
                  integer *ldvt,
                  doublereal *work,
                  integer *lwork,
                  integer *info);

      void cgesvd_(char *jobu,
                  char *jobvt,
                  integer *m,
                  integer *n,
                  complex *a,
                  integer *lda,
                  real *s,
                  complex *u,
                  integer *ldu,
                  complex *vt,
                  integer *ldvt,
                  complex *work,
                  integer *lwork,
                  real *rwork,
                  integer *info);

      void zgesvd_(char *jobu,
                  char *jobvt,
                  integer *m,
                  integer *n,
                  doublecomplex *a,
                  integer *lda,
                  doublereal *s,
                  doublecomplex *u,
                  integer *ldu,
                  doublecomplex *vt,
                  integer *ldvt,
                  doublecomplex *work,
                  integer *lwork,
                  doublereal *rwork,
                  integer *info);

     // cholesky
     void spotrf_(char* uplo,
                 integer* n,
                 real* a,
                 integer* lda,
                 integer* info);

     void cpotrf_(char* uplo,
                 integer* n,
                 complex* a,
                 integer* lda,
                 integer* info);

     void dpotrf_(char* uplo,
                 integer* n,
                 doublereal* a,
                 integer* lda,
                 integer* info);

     void zpotrf_(char* uplo,
                 integer* n,
                 doublecomplex* a,
                 integer* lda,
                 integer* info);

     void spotrs_(char* uplo,
                 integer* n,
                 integer* nrhs,
                 real* a,
                 integer* lda,
                 real* b,
                 integer* ldb,
                 integer* info);

     void cpotrs_(char* uplo,
                 integer* n,
                 integer* nrhs,
                 complex* a,
                 integer* lda,
                 complex* b,
                 integer* ldb,
                 integer* info);

     void dpotrs_(char* uplo,
                 integer* n,
                 integer* nrhs,
                 doublereal* a,
                 integer* lda,
                 doublereal* b,
                 integer* ldb,
                 integer* info);

     void zpotrs_(char* uplo,
                 integer* n,
                 integer* nrhs,
                 doublecomplex* a,
                 integer* lda,
                 doublecomplex* b,
                 integer* ldb,
                 integer* info);

     // LU
      void cgetrf_(integer *m,
                  integer *n,
                  complex *a,
                  integer *lda,
                  integer *ipiv,
                  integer *info);

      void dgetrf_(integer *m,
                  integer *n,
                  doublereal *a,
                  integer *lda,
                  integer *ipiv,
                  integer *info);

      void sgetrf_(integer *m,
                  integer *n,
                  real *a,
                  integer *lda,
                  integer *ipiv,
                  integer *info);

      void zgetrf_(integer *m,
                  integer *n,
                  doublecomplex *a,
                  integer *lda,
                  integer *ipiv,
                  integer *info);

      // inv from LU
      void cgetri_(integer *n,
                  complex *a,
                  integer *lda,
                  integer *ipiv,
                  complex *work,
                  integer *lwork,
                  integer *info);

      void dgetri_(integer *n,
                  doublereal *a,
                  integer *lda,
                  integer *ipiv,
                  doublereal *work,
                  integer *lwork,
                  integer *info);

      void sgetri_(integer *n,
                  real *a,
                  integer *lda,
                  integer *ipiv,
                  real *work,
                  integer *lwork,
                  integer *info);

      void zgetri_(integer *n,
                  doublecomplex *a,
                  integer *lda,
                  integer *ipiv,
                  doublecomplex *work,
                  integer *lwork,
                  integer *info);

     // Perform a series of row interchanges according to a vector of
     // pivot indices
     void claswp_(integer *n,
                 complex *a,
                 integer *lda,
                 integer *k1,
                 integer *k2,
                 integer *ipiv,
                 integer *incx);

     void dlaswp_(integer *n,
                 doublereal *a,
                 integer *lda,
                 integer *k1,
                 integer *k2,
                 integer *ipiv,
                 integer *incx);

     void slaswp_(integer *n,
                 real *a,
                 integer *lda,
                 integer *k1,
                 integer *k2,
                 integer *ipiv,
                 integer *incx);

     void zlaswp_(integer *n,
                 doublecomplex *a,
                 integer *lda,
                 integer *k1,
                 integer *k2,
                 integer *ipiv,
                 integer *incx);


     // Solves a system of linear equations
     //    A * X = B  or  A' * X = B
     // with a general N-by-N matrix A using the LU
     // factorization computed by {C,D,S,Z}GETRF.
     // http://www.netlib.org/lapack/double/dgetrs.f
     // SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

     void cgetrs_(char *trans,
                 integer *n,
                 integer *nrhs,
                 complex *a,
                 integer *lda,
                 integer *ipiv,
                 complex *b,
                 integer *lbd ,
                 integer *info);

     void dgetrs_(char *trans,
                 integer *n,
                 integer *nrhs,
                 doublereal *a,
                 integer *lda,
                 integer *ipiv,
                 doublereal *b,
                 integer *lbd ,
                 integer *info);

     void sgetrs_(char *trans,
                 integer *n,
                 integer *nrhs,
                 real *a,
                 integer *lda,
                 integer *ipiv,
                 real *b,
                 integer *lbd ,
                 integer *info);

     void zgetrs_(char *trans,
                 integer *n,
                 integer *nrhs,
                 doublecomplex *a,
                 integer *lda,
                 integer *ipiv,
                 doublecomplex *b,
                 integer *lbd ,
                 integer *info);

   } // extern "C"
} // namespace fortran

#endif // GLAPACK_HEADERS_H
