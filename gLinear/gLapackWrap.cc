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

char GLAPACKWRAP_CC_VERSION[] = "$Id$";

// glinear headers
#include "gLapackWrap.h"
#include "gRowVector.h"

/* ###################### {S,C,D,Z} POTRF Dispatch ####################### */

// cholesky factorization of a symmetric positive definite matrix
// spotrf dispatch routine (real)
// a [N x N] real input / output
// n rows in a
// info 0 on success
void o__potrf_dispatch(fortran::real* a,
                       fortran::integer n,
                       fortran::integer& info) {
  char uplo = 'L';
  fortran::integer lda = n;

  fortran::spotrf_(&uplo,
                   &n,
                   a,
                   &lda,
                   &info);

}

// cholesky factorization of a symmetric positive definite matrix
// cpotrf dispatch routine (real)
// a [N x N] complex input / output
// n rows in a
// info 0 on success
void o__potrf_dispatch(fortran::complex* a,
                       fortran::integer n,
                       fortran::integer& info) {
  char uplo = 'L';
  fortran::integer lda = n;

  fortran::cpotrf_(&uplo,
                   &n,
                   a,
                   &lda,
                   &info);

}

// cholesky factorization of a symmetric positive definite matrix
// dpotrf dispatch routine (real)
// a [N x N] doublereal input / output
// n rows in a
// info 0 on success
void o__potrf_dispatch(fortran::doublereal* a,
                       fortran::integer n,
                       fortran::integer& info) {
  char uplo = 'L';
  fortran::integer lda = n;

  fortran::dpotrf_(&uplo,
                   &n,
                   a,
                   &lda,
                   &info);

}

// cholesky factorization of a symmetric positive definite matrix
// zpotrf dispatch routine (real)
// a [N x N] doublecomplex input / output
// n rows in a
// info 0 on success
void o__potrf_dispatch(fortran::doublecomplex* a,
                       fortran::integer n,
                       fortran::integer& info) {
  char uplo = 'L';
  fortran::integer lda = n;

  fortran::zpotrf_(&uplo,
                   &n,
                   a,
                   &lda,
                   &info);

}

/* ###################### {S,C,D,Z} GETRF Dispatch ####################### */

// LU factorization
// sgetrf dispatch routine (real)
// a [M x N] real input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output
void o__getrf_dispatch(fortran::real* a,
                       fortran::integer m,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = m;

   fortran::sgetrf_(&m,
        &n,
        a,
        &lda,
        ipiv,
        &info);

}

// LU factorization
// cgetrf dispatch routine (complex)
// a [M x N] complex input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output
void o__getrf_dispatch(fortran::complex* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info) {

   fortran::integer lda = m;

   fortran::cgetrf_(&m,
        &n,
        a,
        &lda,
        ipiv,
        &info);

}

// LU factorization
// dgetrf dispatch routine (doublereal)
// a [M x N] doublereal input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output
void o__getrf_dispatch(fortran::doublereal* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info) {

   fortran::integer lda = m;

   fortran::dgetrf_(&m,
        &n,
        a,
        &lda,
        ipiv,
        &info);

}

// LU factorization
// zgetrf dispatch routine (doublecomplex)
// a [M x N] doublecomplex input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output
void o__getrf_dispatch(fortran::doublecomplex* a,
                       fortran::integer m,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = m;

   fortran::zgetrf_(&m,
        &n,
        a,
        &lda,
        ipiv,
        &info);

}

/* ###################### {S,C,D,Z} POTRS Dispatch ####################### */

// Solve a Linear System using cholesky factorization of a symmetric positive definite matrix
// spotrs dispatch routine (real)
// a [N x N] real input
// n rows in a
// b [Nxnrhs] real rhs input/ output
// nrhs cols in b
// uplo 'L' for lower triangle, 'U' for upper
// info 0 on success

void o__potrs_dispatch(fortran::integer n,
                       fortran::integer nrhs,
                       fortran::real* a,
                       fortran::real* b,
                       fortran::integer& info) {

  char uplo = 'L';
  fortran::integer lda = n;
  fortran::integer ldb = n;

  fortran::spotrs_(&uplo,
                   &n,
                   &nrhs,
                   a,
                   &lda,
                   b,
                   &ldb,
                   &info);

} // o__potrs_dispatch

// Solve a Linear System using cholesky factorization of a symmetric positive definite matrix
// cpotrs dispatch routine (complex)
// a [N x N] complex input
// n rows in a
// b [Nxnrhs] complex rhs input/ output
// nrhs cols in b
// info 0 on success

void o__potrs_dispatch(fortran::integer n,
                       fortran::integer nrhs,
                       fortran::complex* a,
                       fortran::complex* b,
                       fortran::integer& info) {

  char uplo = 'L';
  fortran::integer lda = n;
  fortran::integer ldb = n;

  fortran::cpotrs_(&uplo,
                   &n,
                   &nrhs,
                   a,
                   &lda,
                   b,
                   &ldb,
                   &info);

} // o__potrs_dispatch

// Solve a Linear System using cholesky factorization of a symmetric positive definite matrix
// dpotrs dispatch routine (doublereal)
// a [N x N] doublereal input
// n rows in a
// b [Nxnrhs] doublereal rhs input/ output
// nrhs cols in b
// info 0 on success

void o__potrs_dispatch(fortran::integer n,
                       fortran::integer nrhs,
                       fortran::doublereal* a,
                       fortran::doublereal* b,
                       fortran::integer& info) {

  char uplo = 'L';
  fortran::integer lda = n;
  fortran::integer ldb = n;

  fortran::dpotrs_(&uplo,
                   &n,
                   &nrhs,
                   a,
                   &lda,
                   b,
                   &ldb,
                   &info);

} // o__potrs_dispatch

// Solve a Linear System using cholesky factorization of a symmetric positive definite matrix
// zpotrs dispatch routine (doublecomplex)
// a [N x N] doublecomplex input
// n rows in a
// b [Nxnrhs] doublecomplex rhs input/ output
// nrhs cols in b
// info 0 on success

void o__potrs_dispatch(fortran::integer n,
                       fortran::integer nrhs,
                       fortran::doublecomplex* a,
                       fortran::doublecomplex* b,
                       fortran::integer& info) {

  char uplo = 'L';
  fortran::integer lda = n;
  fortran::integer ldb = n;

  fortran::zpotrs_(&uplo,
                   &n,
                   &nrhs,
                   a,
                   &lda,
                   b,
                   &ldb,
                   &info);

} // o__potrs_dispatch

/* ###################### {S,C,D,Z} GETRS Dispatch ####################### */

// Solve a Linear System using LU factorization
// sgetrs dispatch routine (real)
// a [M x N] real input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output

void o__getrs_dispatch(char trans,
                       fortran::integer n,
                       fortran::integer nrhs,
                       fortran::real* a,
                       fortran::integer* ipiv,
                       fortran::real* b,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer ldb = n;

   fortran::sgetrs_(&trans,
                    &n,
                    &nrhs,
                    a,
                    &lda,
                    ipiv,
                    b,
                    &ldb,
                    &info);

} // o__getrs_dispatch

// Solve a Linear System using LU factorization
// cgetrs dispatch routine (complex)
// a [M x N] real input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output

void o__getrs_dispatch(char trans,
                       fortran::integer n,
                       fortran::integer nrhs,
                       fortran::complex* a,
                       fortran::integer* ipiv,
                       fortran::complex* b,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer ldb = n;

   fortran::cgetrs_(&trans,
                    &n,
                    &nrhs,
                    a,
                    &lda,
                    ipiv,
                    b,
                    &ldb,
                    &info);

} // o__getrs_dispatch

// Solve a Linear System using LU factorization
// dgetrs dispatch routine (doublereal)
// a [M x N] real input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output

void o__getrs_dispatch(char trans,
                       fortran::integer n,
                       fortran::integer nrhs,
                       fortran::doublereal* a,
                       fortran::integer* ipiv,
                       fortran::doublereal* b,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer ldb = n;

   fortran::dgetrs_(&trans,
                    &n,
                    &nrhs,
                    a,
                    &lda,
                    ipiv,
                    b,
                    &ldb,
                    &info);

} // o__getrs_dispatch

// Solve a Linear System using LU factorization
// zgetrs dispatch routine (doublecomplex)
// a [M x N] real input / output
// m rows in a
// n cols in a
// ipiv [1 x min(m,n) ] integer pivot indices output

void o__getrs_dispatch(char trans,
                       fortran::integer n,
                       fortran::integer nrhs,
                       fortran::doublecomplex* a,
                       fortran::integer* ipiv,
                       fortran::doublecomplex* b,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer ldb = n;

   fortran::zgetrs_(&trans,
                    &n,
                    &nrhs,
                    a,
                    &lda,
                    ipiv,
                    b,
                    &ldb,
                    &info);

} // o__getrs_dispatch


/* ###################### {S,C,D,Z} LASWP Dispatch ####################### */
// Apply a permutation
// dlaswp dispatch routine (real)
// a [M x N] doublereal input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
fortran::integer o__laswp_dispatch(fortran::real* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv) {

  //fortran::integer lda = n;
   //fortran::integer lwork = n;
   fortran::integer k1 = 1;
   fortran::integer k2 = m;
   fortran::integer incx = 1;

   fortran::slaswp_(&n,
                    a,
                    &m,
                    &k1,
                    &k2,
                    ipiv,
                    &incx);
   return 1;
}

// Apply a permutation
// dlaswp dispatch routine (complex)
// a [N x N] doublereal input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
fortran::integer o__laswp_dispatch(fortran::complex* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv) {

  //fortran::integer lda = n;
   //fortran::integer lwork = n;
   fortran::integer k1 = 1;
   fortran::integer k2 = m;
   fortran::integer incx = 1;

   fortran::claswp_(&n,
                    a,
                    &m,
                    &k1,
                    &k2,
                    ipiv,
                    &incx);
   return 1;
}

// Apply a permutation
// dlaswp dispatch routine (doublereal)
// a [N x N] doublereal input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
fortran::integer o__laswp_dispatch(fortran::doublereal* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv) {

  //fortran::integer lda = n;
  //fortran::integer lwork = n;
  fortran::integer k1 = 1;
  fortran::integer k2 = m;
  fortran::integer incx = 1;

   fortran::dlaswp_(&n,
                    a,
                    &m,
                    &k1,
                    &k2,
                    ipiv,
                    &incx);
   return 1;
}

// Apply a permutation
// dlaswp dispatch routine (doublecomplex)
// a [N x N] doublereal input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
fortran::integer o__laswp_dispatch(fortran::doublecomplex* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv) {

  //fortran::integer lda = n;
   //fortran::integer lwork = n;
   fortran::integer k1 = 1;
   fortran::integer k2 = m;
   fortran::integer incx = 1;

   fortran::zlaswp_(&n,
                    a,
                    &m,
                    &k1,
                    &k2,
                    ipiv,
                    &incx);
   return 1;
}



/* ###################### {S,C,D,Z} GETRI Dispatch ####################### */

// matrix inversion
// sgetri dispatch routine (real)
// a [N x N] real input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
void o__getri_dispatch(fortran::real* a,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer lwork = n;
   fortran::fortran_array<fortran::real> work(lwork);

   fortran::sgetri_(&n,
        a,
        &lda,
        ipiv,
        work,
        &lwork,
        &info);

}

// matrix inversion
// cgetri dispatch routine (complex)
// a [N x N] complex input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
void o__getri_dispatch(fortran::complex* a,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer lwork = n;
   fortran::fortran_array<fortran::complex> work(lwork);

   fortran::cgetri_(&n,
        a,
        &lda,
        ipiv,
        work,
        &lwork,
        &info);

}

// matrix inversion
// dgetri dispatch routine (double)
// a [N x N] doublereal input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
void o__getri_dispatch(fortran::doublereal* a,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer lwork = n;
   fortran::fortran_array<fortran::doublereal> work(lwork);

   fortran::dgetri_(&n,
        a,
        &lda,
        ipiv,
        work,
        &lwork,
        &info);

}

// matrix inversion
// zgetri dispatch routine (complex)
// a [N x N] doublecomplex input / output
// n order of a (rows/cols)
// ipiv [1 x n ] integer pivot indices output
void o__getri_dispatch(fortran::doublecomplex* a,
                       fortran::integer n,
                       fortran::integer* ipiv,
                       fortran::integer& info) {

   fortran::integer lda = n;
   fortran::integer lwork = n;
   fortran::fortran_array<fortran::doublecomplex> work(lwork);

   fortran::zgetri_(&n,
        a,
        &lda,
        ipiv,
        work,
        &lwork,
        &info);

}


/* ###################### {C,Z} GEEVX Dispatch ####################### */

// cgeevx dispatch routine ( complex<real> )
// a [N x N] complex input data
// n dimension (the N in '[N x N]')
// vr [N x N] complex eigenvector (columns) return
// w [1 x N] complex eigenvalue return
fortran::integer o__geevx_dispatch(fortran::complex* a,
           fortran::integer n,
           fortran::complex* vr,
           fortran::complex* w) {

   char balanc = 'N';   // no balancing
   char jobvl = 'N';     // no left evect's
   char jobvr = 'V';     // do right evect's
   char sense = 'N';     // no reciprocal condition nrs.

   fortran::integer lda = n;  // leading dimension of A
   // fortran::fortran_array<fortran::complex> w(n);  eigenvalue output ( complex array, Nx1)
   fortran::complex vl;      // left-hand Evects (not used here, jobvl == 'N')
   fortran::integer ldvl = 1; // VL leading dim (ignored here, jobvl == 'N')
   // fortran::fortran_array<fortran::complex> vr(n*n);       right-hand Evects.
   fortran::integer ldvr = n;    // VR leading dimension.
   fortran::integer ilo;      // used in balancing
   fortran::integer ihi;       // used in balancing
   fortran::fortran_array<fortran::real> scale(n);      // more balancing detail
   fortran::real abnrm;      // max of sum of abs of elements of each columns
   fortran::fortran_array<fortran::real> rconde(n);    // reciprocal condition nrs. of the Evals
   fortran::fortran_array<fortran::real> rcondv(n);    // reciprocal condition nrs, of the Evects
   fortran::integer lwork = ( (sense == 'N' || sense == 'E' ) ? (2*n) :
               (n*n+2*n) ); // workspace dimension (larger is better)
   fortran::fortran_array<fortran::complex> work(lwork);    // workspace
   fortran::fortran_array<fortran::real> rwork(2*n);      // aux workspace.
   fortran::integer info;      // return value

   fortran::cgeevx_(&balanc,
        &jobvl,
        &jobvr,
        &sense,
        &n,
        a,
        &lda,
        w,
        &vl,
        &ldvl,
        vr,
        &ldvr,
        &ilo,
        &ihi,
        scale,
        &abnrm,
        rconde,
        rcondv,
        work,
        &lwork,
        rwork,
        &info);

   return info;
}

// zgeevx dispatch routine ( complex<double> )
// a [N x N] complex input data
// n dimension (the N in '[N x N]')
// vr [N x N] complex eigenvector (columns) return
// w [1 x N] complex eigenvalue return
fortran::integer o__geevx_dispatch(fortran::doublecomplex* a,
           fortran::integer n,
           fortran::doublecomplex* vr,
           fortran::doublecomplex* w) {

   char balanc = 'N';   // no balancing
   char jobvl = 'N';     // no left evect's
   char jobvr = 'V';     // do right evect's
   char sense = 'N';     // no reciprocal condition nrs.


   fortran::integer lda = n;  // leading dimension of A
   // fortran::fortran_array<fortran::doublecomplex> w(n);  eigenvalue output ( doublecomplex array, Nx1)
   fortran::doublecomplex vl;      // left-hand Evects (not used here, jobvl == 'N')
   fortran::integer ldvl = 1; // VL leading dim (ignored here, jobvl == 'N')
   // fortran::fortran_array<fortran::doublecomplex> vr(n*n);       right-hand Evects.
   fortran::integer ldvr = n;    // VR leading dimension.
   fortran::integer ilo;      // used in balancing
   fortran::integer ihi;       // used in balancing
   fortran::fortran_array<fortran::doublereal> scale(n);      // more balancing detail
   fortran::doublereal abnrm;      // max of sum of abs of elements of each columns
   fortran::fortran_array<fortran::doublereal> rconde(n);    // reciprocal condition nrs. of the Evals
   fortran::fortran_array<fortran::doublereal> rcondv(n);    // reciprocal condition nrs, of the Evects
   fortran::integer lwork = ( (sense == 'N' || sense == 'E' ) ? (2*n) :
               (n*n+2*n) ); // workspace dimension (larger is better)
   fortran::fortran_array<fortran::doublecomplex> work(lwork);    // workspace
   fortran::fortran_array<fortran::doublereal> rwork(2*n);      // aux workspace.
   fortran::integer info;      // return value

   fortran::zgeevx_(&balanc,
        &jobvl,
        &jobvr,
        &sense,
        &n,
        a,
        &lda,
        w,
        &vl,
        &ldvl,
        vr,
        &ldvr,
        &ilo,
        &ihi,
        scale,
        &abnrm,
        rconde,
        rcondv,
        work,
        &lwork,
        rwork,
        &info);

   return info;

}


/* ###################### {S,C,D,Z} GESVD Dispatch ####################### */

// sgesvd dispatch routine ( real )
// a [M x N] real input data
// m a rows (the M in '[M x N]')
// n columns (the N in '[M x N]')
// s real array (dimension min(M,N)
// u [M x M] LH singular matrix
// vt [N x N] transposed RH singular matrix
fortran::integer o__gesvd_dispatch(fortran::real* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::real* u,
           fortran::real* vt) {

   char jobu = 'A';    // return ALL cols. of U in U.
   char jobvt = 'A';  // reutnr ALL cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::real* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::real* s array of diag of S, dimensiom min(M,N)
   // fortran::real* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::real* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = n; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // SGESVD  LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
   fortran::integer lwork = GLINEAR_o__MAX( (3 * GLINEAR_o__MIN(m, n) + GLINEAR_o__MAX(m,n)), (5 * GLINEAR_o__MIN(m,n)) );
   fortran::fortran_array<fortran::real> work(lwork);
   fortran::integer info;

   fortran::sgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        &info);

   return info;
}

// cgesvd dispatch routine ( complex )
// a [M x N] complex input data
// m a rows (the M in '[M x N]')
// n a columns (the N in '[M x N]')
// s real array (dimension min(M,N)
// u [M x M] complex LH singular matrix
// vt [N x N] complex transposed RH singular matrix
fortran::integer o__gesvd_dispatch(fortran::complex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::complex* u,
           fortran::complex* vt) {

   char jobu = 'A';    // return ALL cols. of U in U.
   char jobvt = 'A';  // reutnr ALL cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::real* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::real* s array of diag of S, dimensiom min(M,N)
   // fortran::real* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::real* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = n; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // CGESVD  LWORK >= 2*MIN(M,N)+MAX(M,N)
   fortran::integer lwork = (2 * GLINEAR_o__MIN(m, n)) + GLINEAR_o__MAX(m,n);

   fortran::fortran_array<fortran::complex> work(lwork);
   fortran::fortran_array<fortran::real> rwork(5 * (GLINEAR_o__MAX(m, n)) );
   fortran::integer info;

   fortran::cgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        rwork,
        &info);

   return info;
}

// dgesvd dispatch routine ( doublereal )
// a [M x N] doublereal input data
// m a rows (the M in '[M x N]')
// n columns (the N in '[M x N]')
// s doublereal array (dimension min(M,N)
// u [M x M] doublereal LH singular matrix
// vt [N x N] doublereal transposed RH singular matrix
fortran::integer o__gesvd_dispatch(fortran::doublereal* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublereal* u,
           fortran::doublereal* vt) {

   char jobu = 'A';    // return ALL cols. of U in U.
   char jobvt = 'A';  // reutnr ALL cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::doublereal* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::doublereal* s array of diag of S, dimensiom min(M,N)
   // fortran::doublereal* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::doublereal* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = n; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // DGESVD   LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
   fortran::integer lwork = GLINEAR_o__MAX( (3 * GLINEAR_o__MIN(m, n) + GLINEAR_o__MAX(m,n)), (5 * GLINEAR_o__MIN(m,n)) );
   fortran::fortran_array<fortran::doublereal> work(lwork);
   fortran::integer info;

   fortran::dgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        &info);

   return info;
}

// zgesvd dispatch routine ( doublecomplex )
// a [M x N] doublecomplex input data
// m a rows (the M in '[M x N]')
// n a columns (the N in '[M x N]')
// s real array (dimension min(M,N)
// u [M x M] doublecomplex LH singular matrix
// vt [N x N] doublecomplex transposed RH singular matrix
fortran::integer o__gesvd_dispatch(fortran::doublecomplex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublecomplex* u,
           fortran::doublecomplex* vt) {

   char jobu = 'A';    // return ALL cols. of U in U.
   char jobvt = 'A';  // reutnr ALL cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::doublecomplex* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::doublereal* s array of diag of S, dimensiom min(M,N)
   // fortran::doublecomplexreal* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::doublecomplex* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = n; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // CGESVD  LWORK >= 2*MIN(M,N)+MAX(M,N)
   fortran::integer lwork = (2 * GLINEAR_o__MIN(m, n)) + GLINEAR_o__MAX(m,n);
   fortran::fortran_array<fortran::doublecomplex> work(lwork);
   fortran::fortran_array<fortran::doublereal> rwork(5 * (GLINEAR_o__MAX(m, n)) );
   fortran::integer info;

   fortran::zgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        rwork,
        &info);

   return info;
}


/* ###################### Reduced {S,C,D,Z} GESVD Dispatch ####################### */

// reduced sgesvd dispatch routine ( real )
// a [M x N] real input data
// m a rows (the M in '[M x N]')
// n columns (the N in '[M x N]')
// s real array (dimension min(M,N))
// u - does not return u, returns a dummy
// vt [min(M,N) x N] transposed RH singular matrix
fortran::integer o__reduced_gesvd_dispatch(fortran::real* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::real* u,
           fortran::real* vt) {

   char jobu = 'N';    // return NO cols. of U in U.
   char jobvt = 'S';  // return min(m,n) cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::real* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::real* s array of diag of S, dimensiom min(M,N)
   // fortran::real* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::real* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = m; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // SGESVD  LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
   fortran::integer lwork = GLINEAR_o__MAX( (3 * GLINEAR_o__MIN(m, n) + GLINEAR_o__MAX(m,n)), (5 * GLINEAR_o__MIN(m,n)) );
   fortran::fortran_array<fortran::real> work(lwork);
   fortran::integer info;

   fortran::sgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        &info);

   return info;
}

// reduced cgesvd dispatch r:x
// a [M x N] complex input data
// m a rows (the M in '[M x N]')
// n a columns (the N in '[M x N]')
// s real array (dimension min(M,N))
// u - does not return u, returns a dummy
// vt [min(M,N) x N] transposed RH singular matrix
fortran::integer o__reduced_gesvd_dispatch(fortran::complex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::complex* u,
           fortran::complex* vt) {

   char jobu = 'N';    // return NO cols. of U in U.
   char jobvt = 'S';  // return min(m,n) cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::real* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::real* s array of diag of S, dimensiom min(M,N)
   // fortran::real* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::real* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = m; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // CGESVD  LWORK >= 2*MIN(M,N)+MAX(M,N)
   fortran::integer lwork = (2 * GLINEAR_o__MIN(m, n)) + GLINEAR_o__MAX(m,n);

   fortran::fortran_array<fortran::complex> work(lwork);
   fortran::fortran_array<fortran::real> rwork(5 * (GLINEAR_o__MAX(m, n)) );
   fortran::integer info;

   fortran::cgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        rwork,
        &info);

   return info;
}

// reduced dgesvd dispatch routine ( doublereal )
// a [M x N] doublereal input data
// m a rows (the M in '[M x N]')
// n columns (the N in '[M x N]')
// s doublereal array (dimension min(M,N))
// u - does not return u, returns a dummy
// vt [min(M,N) x N] transposed RH singular matrix
fortran::integer o__reduced_gesvd_dispatch(fortran::doublereal* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublereal* u,
           fortran::doublereal* vt) {

   char jobu = 'N';    // return NO cols. of U in U.
   char jobvt = 'S';  // return min(m,n) cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::doublereal* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::doublereal* s array of diag of S, dimensiom min(M,N)
   // fortran::doublereal* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::doublereal* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = m; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // DGESVD   LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
   fortran::integer lwork = GLINEAR_o__MAX( (3 * GLINEAR_o__MIN(m, n) + GLINEAR_o__MAX(m,n)), (5 * GLINEAR_o__MIN(m,n)) );
   fortran::fortran_array<fortran::doublereal> work(lwork);
   fortran::integer info;

   fortran::dgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        &info);

   return info;
}

// reduced zgesvd dispatch routine ( doublecomplex )
// a [M x N] doublecomplex input data
// m a rows (the M in '[M x N]')
// n a columns (the N in '[M x N]')
// s real array (dimension min(M,N))
// u - does not return u, returns a dummy
// vt [min(M,N) x N] transposed RH singular matrix
fortran::integer o__reduced_gesvd_dispatch(fortran::doublecomplex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublecomplex* u,
           fortran::doublecomplex* vt) {

   char jobu = 'N';    // return NO cols. of U in U.
   char jobvt = 'S';  // return min(m,n) cols. of V' in VT.
   // fortran::integer m   number of rows in A
   // fortran::integer n    number or cols in A
   // fortran::doublecomplex* a [M x N] input matrix
   fortran::integer lda = m; // leading dimension of A
   // fortran::doublereal* s array of diag of S, dimensiom min(M,N)
   // fortran::doublecomplexreal* u [ LDU x M ] LH singlular matrix
   fortran::integer ldu = m; // leading dimension of U = M
   // fortran::doublecomplex* vt [ LDVT x N ] RH singluar matrix
   fortran::integer ldvt = m; // leading dimension of VT = N
   // fortran::integer lwork = o__MAX( (3 * o__MAX(m, n) + o__MAX(m,n)), (5 * o__MIN(m,n) - 4) );
   // CGESVD  LWORK >= 2*MIN(M,N)+MAX(M,N)
   fortran::integer lwork = (2 * GLINEAR_o__MIN(m, n)) + GLINEAR_o__MAX(m,n);
   fortran::fortran_array<fortran::doublecomplex> work(lwork);
   fortran::fortran_array<fortran::doublereal> rwork(5 * (GLINEAR_o__MAX(m, n)) );
   fortran::integer info;

   fortran::zgesvd_(&jobu,
        &jobvt,
        &m,
        &n,
        a,
        &lda,
        s,
        u,
        &ldu,
        vt,
        &ldvt,
        work,
        &lwork,
        rwork,
        &info);

   return info;
}

/*
  Support functions for pivot to permutation convertion
 */
void showipiv(const fortran::fortran_array<fortran::integer>& ipiv,
              int sz) {
  for (fortran::integer i = 0; i < sz; i++) {
    std::cout << ipiv[i] << " ";
  } // for i
  std::cout << std::endl;
} // showipiv

/**
 * transforms lapack ipiv to row permutation vector
 * both arrays MUST be presized
 */
void ipiv2perm(const fortran::fortran_array<fortran::integer>& ipiv,
               gLinear::gRowVector<int>& perm) {
  perm.set(0,1);
  fortran::integer sz = perm.size();
  for (fortran::integer i = 0;  i < sz; i++) {
    std::swap(perm[ ipiv[i]-1 ], perm[i]);
  } // for
} // ipiv2perm

/**
 * transforms row permutation vector to lapack ipiv
 * both arrays must be presized
 */
void perm2ipiv(const gLinear::gRowVector<int>& perm,
               fortran::fortran_array<fortran::integer>& ipiv) {
  fortran::integer sz = perm.size();
  while (sz--) {
    ipiv[sz] = sz;
  } // while

  sz = perm.size();
  for (fortran::integer i = 0;  i < sz; i++) {
    ipiv[i] = ipiv[ perm[i] ];
  } // for

  while(sz--) {
    ipiv[sz]++;
  } // for

} // perm2ipiv

