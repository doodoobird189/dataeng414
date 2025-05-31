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

// gLinear headers
#include "gCppTraits.h"

// standard headers
#include <complex>  // complex, real, imag
#include <iostream>  // cout, endl


GLINEAR_NAMESPACE_BEGIN

/*
 * Author's Disclaimer:
 *
 * The following code is not meant to serve as an example of good programming practice.
 * It should in fact not serve as an example of anything.
 * So, if you read through this thinking that it is the worst bowl of spaghetti you've
 * ever seen, just remember: you weren't the first to think so.
 */


template<typename TVect, typename TVal>
inline gEigenData<TVect, TVal>::gEigenData(const gMatrix<T_EVector, ROW_DENSE>& newvect,
             const gMatrix<T_EValue, ROW_DENSE>& newval) :
  eVectors(newvect),
  eValues(newval) {
} // constructor

template<typename TVect, typename TVal>
template<typename T_srcVec, typename T_srcVal>
inline gEigenData<TVect, TVal>&
gEigenData<TVect, TVal>::operator=(const gEigenData<T_srcVec, T_srcVal>& src) {
  eVectors = src.vectors();
  eValues = src.values();
  return *this;
} // operator=

template<typename TVect, typename TVal>
inline gMatrix<typename gEigenData<TVect, TVal>::T_EVector, ROW_DENSE>&
gEigenData<TVect, TVal>::vectors(void) {
  return eVectors;
} // vectors

template<typename TVect, typename TVal>
inline const gMatrix<typename gEigenData<TVect, TVal>::T_EVector, ROW_DENSE>&
gEigenData<TVect, TVal>::vectors(void) const {
  return eVectors;
} // vectors

template<typename TVect, typename TVal>
inline gMatrix<typename gEigenData<TVect, TVal>::T_EValue, ROW_DENSE>&
gEigenData<TVect, TVal>::values(void) {
  return eValues;
} // values

template<typename TVect, typename TVal>
inline const gMatrix<typename gEigenData<TVect, TVal>::T_EValue, ROW_DENSE>&
gEigenData<TVect, TVal>::values(void) const {
  return eValues;
} // values



/* #################### Support Functions ############  */

#define GLINEAR_o__MAX(a,b) ( (a > b) ? (a) : (b) )
#define GLINEAR_o__MIN(a,b) ( (a < b) ? (a) : (b) )

// these work like UNIX copy command (cp): source first, then destination

// complex<T>, gFortranTraits<T>::T_Complex
template<typename T>
inline void o__copy_c_to_f(const std::complex<T>& src,
                           typename gFortranTraits<T>::T_Complex& dst) {
  dst.r = std::real(src);
  dst.i = std::imag(src);
} // o__copy_c_to_f

// T, gFortranTraits<T>::T_Complex
template<typename T>
inline void o__copy_c_to_f(const T& src,
                           typename gFortranTraits<T>::T_Complex& dst) {
  dst.r = src;
  dst.i = 0;
} // o__copy_c_to_f

// T, gFortranTraits<T>::T_Fortran
template<typename T>
inline void o__copy_c_to_f(const T& src,
                           typename gFortranTraits<T>::T_Fortran& dst) {
  dst = src;
} // o__copy_c_to_f

// gFortranTraits<T>::T_Complex, T
template<typename T>
inline void o__copy_f_to_c(const typename gFortranTraits<T>::T_Complex& src,
                           std::complex<T>& dst) {
  dst = std::complex<T>(src.r, src.i);
} // o__copy_f_to_c

// gFortranTraits<T>::T_Fortran, T
template<typename T>
inline void o__copy_f_to_c(const typename gFortranTraits<T>::T_Fortran& src,
                           T& dst) {
  dst = src;
} // o__copy_f_to_c

// gFortranTraits<T>::T_Fortran -> T
template<typename T>
inline T o__conv_f_to_c(const typename gFortranTraits<T>::T_Fortran& src) {
  return src;
} // o__conv_f_to_c

// gFortranTraits<T>::T_Complex -> std::complex<T>
template<typename T>
inline std::complex<T> o__conv_f_to_c(const typename gFortranTraits<T>::T_Complex& src) {
  return std::complex<T>(src.r, src.i);
} // o__conv_f_to_c

// copy C++ gColVector to corresponding Fortran array. copy T1 -> T2.
// Both structures must be presized.
template<typename T1, typename T2>
inline void o__copy_c_to_f(const gColVector<T1>& src,
                           fortran::fortran_array<T2>& dst) {
  typename gColVector<T1>::T_Index sz = src.size();
  for (T_Index indx = 0; indx < sz; indx++) {
      o__copy_c_to_f( src[indx], dst[indx] );
  } // for
} // o__copy_c_to_f

// copy Fortran array to corresponding C++ gColVector. copy T2 -> T1.
// Both structures must be presized.
template<typename T1, typename T2>
inline void o__copy_f_to_c(const fortran::fortran_array<T2>& src,
                           gColVector<T1>& dst) {
  typename gColVector<T1>::T_Index sz = dst.size();
  for (T_Index indx = 0; indx < sz; indx++) {
      o__copy_f_to_c( src[indx], dst[indx] );
  } // for
} // o__copy_f_to_c

// copy C++ gRowVector to corresponding Fortran array. copy T1 -> T2.
// Both structures must be presized.
template<typename T1, typename T2>
inline void o__copy_c_to_f(const gRowVector<T1>& src,
                           fortran::fortran_array<T2>& dst) {
  typename gRowVector<T1>::T_Index sz = src.size();
  for (T_Index indx = 0; indx < sz; indx++) {
      o__copy_c_to_f( src[indx], dst[indx] );
  } // for
} // o__copy_c_to_f

// copy Fortran array to corresponding C++ gRowVector. copy T2 -> T1.
// Both structures must be presized.
template<typename T1, typename T2>
inline void o__copy_f_to_c(const fortran::fortran_array<T2>& src,
                           gRowVector<T1>& dst) {
  typename gRowVector<T1>::T_Index sz = dst.size();
  for (T_Index indx = 0; indx < sz; indx++) {
      o__copy_f_to_c( src[indx], dst[indx] );
  } // for
} // o__copy_f_to_c

// copy C++ Matrix to corresponding Fortran array. copy T1 -> T2.
// Both structures must be presized.
template<typename T1, typename T_Storage, typename T2>
inline void o__copy_c_to_f(const gMatrix<T1, T_Storage>& src,
                           fortran::fortran_array<T2>& dst) {
  typename gMatrix<T1, T_Storage>::T_Index nr = src.rows();
  typename gMatrix<T1, T_Storage>::T_Index nc = src.cols();
  typedef typename gMatrix<T1, T_Storage>::T_Index T_Index;
  T_Index indx = 0;
  for (T_Index cindx = 0; cindx < nc; cindx++) {
    for (T_Index rindx = 0; rindx < nr; rindx++) {
      o__copy_c_to_f( src(rindx, cindx), dst[indx++] );
    } // for
  } // for
} // o__copy_c_to_f

// copy Fortran array to corresponding C++ gMatrix. copy T2 -> T1.
// Both structures must be presized.
template<typename T1, typename T_Storage, typename T2>
inline void o__copy_f_to_c(const fortran::fortran_array<T2>& src,
                           gMatrix<T1, T_Storage>& dst) {
  typename gMatrix<T1, T_Storage>::T_Index nr = dst.rows();
  typename gMatrix<T1, T_Storage>::T_Index nc = dst.cols();
  typedef typename gMatrix<T1, T_Storage>::T_Index T_Index;
  T_Index indx = 0;
  for (T_Index cindx = 0; cindx < nc; cindx++) {
    for (T_Index rindx = 0; rindx < nr; rindx++) {
      o__copy_f_to_c( src[indx++], dst(rindx, cindx) );
    } // for
  } // for
} // o__copy_f_to_c

// copy Fortran array to corresponding C++ gMatrix. copy T2 -> T1.
// Both structures must be precised
template<typename T1, typename T_Storage, typename T2>
inline void o__copylower_f_to_c(const fortran::fortran_array<T2>& src,
                           gMatrix<T1, T_Storage>& dst) {
  typename gMatrix<T1, T_Storage>::T_Index nr = dst.rows();
  typename gMatrix<T1, T_Storage>::T_Index nc  = dst.cols();
  typedef typename gMatrix<T1, T_Storage>::T_Index T_Index;
  T_Index indx = 0;
  // zero the upper triangle while also moving the src index over the
  // upper triangle rows
  for (T_Index cindx = 0; cindx < nc; cindx++) {
    for (T_Index rindx = 0; rindx < cindx; rindx++) {
      dst(rindx, cindx) = 0; indx++;
    } // for
    for (T_Index rindx = cindx; rindx < nr; rindx++) {
      o__copy_f_to_c( src[indx++], dst(rindx, cindx) );
    } // for
  } // for
} // o__copylower_f_to_c

/* #################### Norms and Condition Numbers ###########  */

// P-norm: only 2-norm (default) currently implemented
template<typename T, typename T_Storage>
inline T matrixNorm(const gMatrix<T, T_Storage>& A,
        int P) {
  if (P == 2) {
    return matrixNorm2(A);
  } // if
#ifdef CHK_IDX
  gErrorInvalidArgument("", "T matrixNorm(const gMatrix<T, T_Storage>& A, int P)", "Only 2-norm supported");
#endif // CHK
  return 0;
} // matrixNorm

// 2-norm: maximum singular value
template<typename T, typename T_Storage>
T matrixNorm2(const gMatrix<T, T_Storage>& A) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef typename gMatrix<T, ROW_DENSE>::T_Index T_Index;
  T_Index M = A.rows();
  T_Index N = A.cols();
  gMatrix<T, ROW_DENSE> ev(N,N);
  gMatrix<T, ROW_DENSE> U(M, M);
  gMatrix<T_Real, ROW_DENSE> S(M, N);
  gMatrix<T, ROW_DENSE> VT(N, N);
  svd(A, U, S, VT);
  T_Index minMN = ( ( M < N ) ? M : N );
  T maxval = S(0,0);
  for (T_Index i = 1; i < minMN; i++) {
    if (S(i,i) > maxval) {
      maxval = S(i,i);
    } // if
  } // for
  return maxval;
} // matrixNorm2

// return P-norm condition number: norm(A, P) * norm(inv(A), P)
// currently only supports norm-2 (default)
template<typename T, typename T_Storage>
T cond(const gMatrix<T, T_Storage>& A,
       int P) {
  if (P == 2) {
    return cond2(A);
  } // if
#ifdef CHK_IDX
  gErrorInvalidArgument("", "T cond(const gMatrix<T, T_Storage>& A, int P)", "Only 2-norm supported");
#endif // CHK_IDX
  return 0;
} // cond

// return 2-norm condition number: max(svd(A)) / min(svd(A))
template<typename T, typename T_Storage>
T cond2(const gMatrix<T, T_Storage>& A) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef typename gMatrix<T, ROW_DENSE>::T_Index T_Index;
  T_Index M = A.rows();
  T_Index N = A.cols();
  gMatrix<T, ROW_DENSE> U(M, M);
  gMatrix<T_Real, ROW_DENSE> S(M, N);
  gMatrix<T, ROW_DENSE> VT(N, N);
  svd(A, U, S, VT);
  T_Index minMN = ( ( M < N ) ? M : N );
  T maxval = S(0,0);
  T minval = S(0,0);
  for (T_Index i = 1; i < minMN; i++) {
    if (S(i,i) > maxval) {
      maxval = S(i,i);
    } // if
    if (S(i,i) < minval) {
      minval = S(i,i);
    } // if
  } // for
  return maxval / minval;
} // cond2

/* #################### Matrix Inversion ############  */

// user-callable inv.
// implementation for general matrices.
template<typename T, typename T_Storage>
inline typename gMatrix<T, T_Storage>::T_Dense
inv(const gMatrix<T, T_Storage>& A,
    int& fail) {
  return o__gentype_inv(A, fail);
} // inv

// user-callable inv.
// implementation for general matrices.
template<typename T, typename T_Storage>
inline typename gMatrix<T, T_Storage>::T_Dense
inv(const gMatrix<T, T_Storage>& A,
    T& det,
    int& fail) {
  return o__gentype_inv(A, det, fail);
} // inv

// user-callable inv.
// implementation for general matrices.
template<typename T, typename T_Storage1, typename T_Storage2>
inline void
inv(const gMatrix<T, T_Storage1>& A,
    gMatrix<T, T_Storage2>& invA,
    int& fail) {
  o__gentype_inv(A, invA, fail);
} // inv

// user-callable inv.
// implementation for general matrices.
template<typename T, typename T_Storage1, typename T_Storage2>
inline void
inv(const gMatrix<T, T_Storage1>& A,
    gMatrix<T, T_Storage2>& invA,
    T& det,
    int& fail) {
  o__gentype_inv(A, invA, det, fail);
} // inv

// helper function that calls gentype_inv returning
// with an arg. matrix
template<typename T, typename T_Storage>
inline typename gMatrix<T, T_Storage>::T_Dense
o__gentype_inv(const gMatrix<T, T_Storage>& A,
               int& fail) {
  typename gMatrix<T, T_Storage>::T_Dense invA;
  o__gentype_inv(A, invA, fail);
  return invA;
} // o__gentype_inv

// helper function that calls gentype_inv returning
// with an arg. matrix
template<typename T, typename T_Storage>
inline typename gMatrix<T, T_Storage>::T_Dense
o__gentype_inv(const gMatrix<T, T_Storage>& A,
               T& det,
               int& fail) {
  typename gMatrix<T, T_Storage>::T_Dense invA;
  o__gentype_inv(A, invA, det, fail);
  return invA;
} // o__gentype_inv

// calc matrix inverse
// helper function that calls the appropriate Lapack functions
// first calls getrf to LU factorize A, then getri to invert from
// LU matrix. copies A as the Lapack trashes it.
template<typename T, typename T_Storage1, typename T_Storage2>
void o__gentype_inv(const gMatrix<T, T_Storage1>& A,
                    gMatrix<T, T_Storage2>& invA,
                    int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
     "Matrix inversion",
     "Matrix not square",
     A.rows(),
     A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer N = A.rows();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
  fortran::fortran_array<fortran::integer> ipiv(N);
  fortran::integer info;

  invA.resize(N, N);
  o__copy_c_to_f(A, a);
  // lu
  o__getrf_dispatch(a, N, N, ipiv, info);
  if (!info) {
    // inv
    o__getri_dispatch(a, N, ipiv, info);
  } // if not fail
  if (!info) {
    o__copy_f_to_c(a, invA);
  } // if not fail
  fail = info;
} // o__gentype_inv

// calc matrix inverse and determinant
// helper function that calls the appropriate Lapack functions
// first calls getrf to LU factorize A, then getri to invert from
// LU matrix. copies A as the Lapack trashes it.
template<typename T, typename T_Storage1, typename T_Storage2>
void o__gentype_inv(const gMatrix<T, T_Storage1>& A,
                    gMatrix<T, T_Storage2>& invA,
                    T& det,
                    int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
     "Matrix inversion",
     "Matrix not square",
     A.rows(),
     A.cols() );
  } // if
#endif // CHK_IDX
  // typedef typename gMatrix<T, T_Storage1>::T_Index T_Index;
  typedef typename gFortranTraits<T>::T_Real T_Real;

  fortran::integer N = A.rows();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
  fortran::fortran_array<fortran::integer> ipiv(N);
  fortran::integer info;

  invA.resize(N, N);
  o__copy_c_to_f(A, a);
  // lu
  o__getrf_dispatch(a, N, N, ipiv, info);
  if (!info) {
    // calc. the det by multiplying the diagonal elements of a, i.e. diagonal of U
    det = o__conv_f_to_c<T_Real>(a[0]);
    fortran::integer indx = N;
    // true => det(P) = -1, false => det(P) = 1
    bool sign = (ipiv[0] != 1);
    while (--indx) {
      det *= o__conv_f_to_c<T_Real>( a[ indx * (N + 1) ] );
      if (ipiv[indx] != indx + 1) {
        sign = !sign;
      } // if
    } // while
    if (sign) {
      det *= static_cast<T>(-1);
    } // if
    // 1 / det(A) = det(inv(A))
    det = 1 / det;
    // inv
    o__getri_dispatch(a, N, ipiv, info);
  } // if not fail
  if (!info) {
    o__copy_f_to_c(a, invA);
  } // if not fail
  fail = info;
} // o__gentype_inv

/* #################### Determinant ############  */

// determinant of general matrix
template<typename T, typename T_Storage>
inline T det(const gMatrix<T, T_Storage>& A,
             int& fail) {
  return o__gentype_det(A, fail);
} // det

// helper function that calls Lapack dispatch and calcs det
// with an arg. matrix
template<typename T, typename T_Storage>
inline T
o__gentype_det(const gMatrix<T, T_Storage>& A,
               int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
     "Matrix determinant",
     "Matrix not square",
     A.rows(),
     A.cols() );
  } // if
#endif // CHK_IDX

  typedef typename gFortranTraits<T>::T_Real T_Real;

  T det;
  fortran::integer N = A.rows();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
  fortran::fortran_array<fortran::integer> ipiv(N);
  fortran::integer info;

  o__copy_c_to_f(A, a);
  // lu
  o__getrf_dispatch(a, N, N, ipiv, info);
  if (!info) {
    // calc. the det by multiplying the diagonal of a, i.e. diagonal of U
    det = o__conv_f_to_c<T_Real>(a[0]);
    fortran::integer indx = N;
    // true => det(P) = -1, false => det(P) = 1
    bool sign = (ipiv[0] != 1);
    while (--indx) {
      det *= o__conv_f_to_c<T_Real>( a[ indx * (N + 1) ] );
      if (ipiv[indx] != indx + 1) {
        sign = !sign;
      } // if
    } // while
    if (sign) {
      det *= static_cast<T>(-1);
    } // if
  } // if not fail
  fail = info;
  return det;
} // o__gentype_det

/* #################### LU Factorization ############  */

// lu of general matrix
// template<typename T, typename T_Storage1, typename T_Storage2>
// inline void lu(const gMatrix<T, T_Storage1>& A,
//                gMatrix<T, T_Storage2>& LU) {
//   o__gentype_lu(A, LU);
// } // lu

// // helper function that calls Lapack dispatch and calcs LU
// // factorization
// template<typename T, typename T_Storage1, typename T_Storage2>
// inline void
// o__gentype_lu(const gMatrix<T, T_Storage1>& A,
//               gMatrix<T, T_Storage2>& LU) {
// #ifdef CHK_IDX
//   if ( A.rows() != A.cols() ) {
//     gErrorShape( "",
//      "Matrix LU",
//      "Matrix not square",
//      A.rows(),
//      A.cols() );
//   } // if
// #endif // CHK_IDX

//   fortran::integer N = A.rows();
//   fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
//   fortran::fortran_array<fortran::integer> ipiv(N);

//   LU.resize(N, N);
//   o__copy_c_to_f(A, a);
//   // lu
//   o__getrf_dispatch(a, N, N, ipiv);
//   o__copy_f_to_c(a, LU);
// } // o__gentype_lu

/**
 * Find the lower triangle cholesky factor of a symmetrical A
 */
template<typename T, typename T_Storage1, typename T_Storage2>
inline void cholesky_factor(const gMatrix<T, T_Storage1>& A,
                            gMatrix<T, T_Storage2>& L,
                            int& fail) {

  o__gentype_cholesky_factor(A, L, fail);

} // cholesky_factor

/**
 * Find the lower triangle cholesky factor of a symmetrical A.
 * The result is kept in the provided Fortran structure to be
 * re-used by the Cholesky solvers that accept a precomputed
 * Cholesky factor.
 */
template<typename T, typename T_Storage1>
inline void cholesky_factor(const gMatrix<T, T_Storage1>& A,
                            fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                            int& fail) {

  o__gentype_cholesky_factor(A, L, fail);

} // cholesky_factor

// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1, typename T_Storage2>
inline void o__gentype_cholesky_factor(const gMatrix<T, T_Storage1>& A,
                                       gMatrix<T, T_Storage2>& L,
                                       int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
                 "Matrix choleski",
                 "Matrix not square",
                 A.rows(),
                 A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer info;
  fortran::integer N = A.rows();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> chol(N*N);
  L.resize(N,N);

  o__copy_c_to_f(A, chol);
  o__potrf_dispatch(chol, N, info);
  if (!info) {
    o__copylower_f_to_c(chol, L);
  } // if not fail
  fail = info;
} // o__gentype_cholesky_factor

// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1>
inline void o__gentype_cholesky_factor(const gMatrix<T, T_Storage1>& A,
                                       fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                                       int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
                 "Matrix choleski",
                 "Matrix not square",
                 A.rows(),
                 A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer info;
  fortran::integer N = A.rows();

  o__copy_c_to_f(A, L);
  o__potrf_dispatch(L, N, info);
  fail = info;
} // o__gentype_cholesky_factor


// user-callable lu factorization with permutation information
template<typename T, typename T_Storage1, typename T_Storage2>
inline void lu_factor(const gMatrix<T, T_Storage1>& A,
                      gMatrix<T, T_Storage2>& LU,
                      gRowVector<int>& P,
                      int& fail) {

  o__gentype_lu_factor(A, LU, P, fail);

} // lu_factor

// helper function that calls Lapack dispatch and calcs LU
// factorization
template<typename T, typename T_Storage1, typename T_Storage2>
inline void
o__gentype_lu_factor(const gMatrix<T, T_Storage1>& A,
                     gMatrix<T, T_Storage2>& LU,
                     gRowVector<int>& P,
                     int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
     "Matrix LU",
     "Matrix not square",
     A.rows(),
     A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer info;
  fortran::integer N = A.rows();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
  fortran::fortran_array<fortran::integer> ipiv(N);
  LU.resize(N, N);
  P.resize(N);

  o__copy_c_to_f(A, a);

  // lu
  o__getrf_dispatch(a, N, N, ipiv, info);
  if (!info) {
    o__copy_f_to_c(a, LU);
    ipiv2perm(ipiv, P);
  } // if not fail
  fail = info;
} // o__gentype_lu_factor

// apply permutation to matrix
template<typename T, typename T_Storage1>
inline void lu_permute(gMatrix<T, T_Storage1>& A,
                       const gRowVector<int>& P) {

  fortran::integer M = A.rows();
  fortran::integer N = A.cols();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(M*N);
  fortran::fortran_array<fortran::integer> ipiv(M);

  o__copy_c_to_f(A, a);
  perm2ipiv(P, ipiv);
  o__laswp_dispatch(a, M, N, ipiv);
  o__copy_f_to_c(a, A);

} // lu_permute


/* #################### Linear Systems ############  */

/**
 * Solve a Single Linear System with symmetrical A
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void cholesky_solve(const gMatrix<T, T_Storage1>& A,
                    const gColVector<T, T_Storage2>& b,
                    gColVector<T, T_Storage3>& x,
                    int& fail) {

  const gMatrix<T> bmat(b);
  gMatrix<T> xmat(x);
  o__gentype_cholesky_solve(A, bmat, xmat, fail);

} // cholesky_solve

/**
 * Solves Multiple Linear Systems with symmetrical A
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void cholesky_solve(const gMatrix<T, T_Storage1>& A,
                    const gMatrix<T, T_Storage2>& B,
                    gMatrix<T, T_Storage3>& X,
                    int& fail){

  o__gentype_cholesky_solve(A, B, X, fail);

} // cholesky_solve

/**
 * Solve a Single Linear System with precomputed Cholesky factor for a symmetrical A
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void cholesky_solve(fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                    const gColVector<T, T_Storage1>& b,
                    gColVector<T, T_Storage2>& x,
                    int& fail) {

  const gMatrix<T> bmat(b);
  gMatrix<T> xmat(x);
  o__gentype_cholesky_solve(L, bmat, xmat, fail);

} // cholesky_solve

/**
 * Solve Multiple Linear System with precomputed Cholesky factor for a symmetrical A
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void cholesky_solve(fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                    const gMatrix<T, T_Storage1>& B,
                    gMatrix<T, T_Storage2>& X,
                    int& fail) {

  o__gentype_cholesky_solve(L, B, X, fail);

} // cholesky_solve

/**
 * Solve a Single Linear System
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& A,
              const gColVector<T, T_Storage2>& b,
              gColVector<T, T_Storage3>& x,
              int& fail,
              const char trans) {

  const gMatrix<T> bmat(b);
  gMatrix<T> xmat(x);
  o__gentype_lu_solve(A, bmat, xmat, trans, fail);

} // lu_solve

/**
 * Solve multiple Linear Systems
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& A,
              const gMatrix<T, T_Storage2>& b,
              gMatrix<T, T_Storage3>& x,
              int& fail,
              const char trans) {

  o__gentype_lu_solve(A, b, x, trans, fail);

} // lu_solve

/**
 * Solve multiple Linear Systems with provided LU factorization
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& LU,
              const gRowVector<int>& P,
              const gMatrix<T, T_Storage2>& b,
              gMatrix<T, T_Storage3>& x,
              int& fail,
              const char trans) {

  o__gentype_lu_solve(LU, P, b, x, trans, fail);

} // lu_solve

/**
 * Solve a Linear System with provided LU factorization
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& LU,
              const gRowVector<int>& P,
              const gColVector<T, T_Storage2>& b,
              gColVector<T, T_Storage3>& x,
              int& fail,
              const char trans) {

  const gMatrix<T> bmat(b);
  gMatrix<T> xmat(x);
  o__gentype_lu_solve(LU, P, bmat, xmat, trans, fail);

} // lu_solve


// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
inline void o__gentype_cholesky_solve(const gMatrix<T, T_Storage1>& A,
                                      const gMatrix<T, T_Storage2>& b,
                                      gMatrix<T, T_Storage3>& x,
                                      int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
                 "Matrix choleski",
                 "Matrix not square",
                 A.rows(),
                 A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer info;
  fortran::integer N = A.rows();
  fortran::integer nrhs = b.cols();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> chol(N*N);
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> rhs(N*nrhs);

  x.resize(N, nrhs);
  o__copy_c_to_f(A, chol);
  o__copy_c_to_f(b, rhs);
  // chol
  o__potrf_dispatch(chol, N, info);
  // solve
  if (!info) {
    o__potrs_dispatch(N, nrhs, chol, rhs,info);
  } // if not fail
  if (!info) {
    o__copy_f_to_c(rhs, x);
  } // if not fail
  fail = info; //o__copy_f_to_c(info, fail);
} // o__gentype_cholesky_solve

// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1, typename T_Storage2>
inline void o__gentype_cholesky_solve(fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                                      const gMatrix<T, T_Storage1>& b,
                                      gMatrix<T, T_Storage2>& x,
                                      int& fail) {

  fortran::integer info;
  fortran::integer N = b.rows();
  fortran::integer nrhs = b.cols();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> rhs(N*nrhs);

  x.resize(N, nrhs);
  o__copy_c_to_f(b, rhs);
  // solve
  o__potrs_dispatch(N, nrhs, L, rhs,info);

  if (!info) {
    o__copy_f_to_c(rhs, x);
  } // if not fail
  fail = info;
} // o__gentype_cholesky_solve


// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
inline void o__gentype_lu_solve(const gMatrix<T, T_Storage1>& A,
                                const gMatrix<T, T_Storage2>& b,
                                gMatrix<T, T_Storage3>& x,
                                const char trans,
                                int& fail) {
#ifdef CHK_IDX
  if ( A.rows() != A.cols() ) {
    gErrorShape( "",
                 "Matrix LU",
                 "Matrix not square",
                 A.rows(),
                 A.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer N = A.rows();
  fortran::integer nrhs = b.cols();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> LU(N*N);
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> rhs(N*nrhs);
  fortran::fortran_array<fortran::integer> ipiv(N);
  fortran::integer info;

  x.resize(N, nrhs);
  o__copy_c_to_f(A, LU);
  o__copy_c_to_f(b, rhs);
  // lu
  o__getrf_dispatch(LU, N, N, ipiv, info);
  // solve
  if (!info) {
    o__getrs_dispatch(trans, N, nrhs, LU, ipiv, rhs, info);
  } // if not fail
  if (!info) {
    o__copy_f_to_c(rhs, x);
  } // if not fail
  fail = info;
} // o__gentype_lu_solve


// helper function that calls Lapack dispatch
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
inline void o__gentype_lu_solve(const gMatrix<T, T_Storage1>& LU,
                                const gRowVector<int>& P,
                                const gMatrix<T, T_Storage2>& b,
                                gMatrix<T, T_Storage3>& x,
                                const char trans,
                                int& fail) {
#ifdef CHK_IDX
  if ( LU.rows() != LU.cols() ) {
    gErrorShape( "",
                 "Matrix LU",
                 "Matrix not square",
                 LU.rows(),
                 LU.cols() );
  } // if
#endif // CHK_IDX

  fortran::integer N = LU.rows();
  fortran::integer nrhs = b.cols();
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> a(N*N);
  fortran::fortran_array<typename gFortranTraits<T>::T_Fortran> rhs(N*nrhs);
  fortran::fortran_array<fortran::integer> ipiv(N);
  fortran::integer info;

  x.resize(N, nrhs);
  o__copy_c_to_f(LU, a);
  o__copy_c_to_f(b, rhs);
  perm2ipiv(P, ipiv);
  // solve
  o__getrs_dispatch(trans, N, nrhs, a, ipiv, rhs, info);
  if (!info) {
    o__copy_f_to_c(rhs, x);
  } // if not fail
  fail = info;
} // o__gentype_lu_solve


/* #################### Eigenanalysis ############  */

// user-callable: EigenData eig ( T )
template<typename T, typename T_Storage>
gEigenData< std::complex<T>, std::complex<T> >
eig(const gMatrix<T, T_Storage>& data) {
  return o__gentype_eig(data);
} // eig

// user-callable: void eig ( T )
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void eig(const gMatrix<T, T_Storage1>& data,
         gMatrix<std::complex<T>, T_Storage2>& evectors,
         gMatrix<std::complex<T>, T_Storage3>& evalues) {
  o__gentype_eig(data, evectors, evalues);
} // eig

// user-callable: gEigenData eig( std::complex<T> )
template<typename T, typename T_Storage>
gEigenData< std::complex<T>, std::complex<T> >
eig(const gMatrix<std::complex<T>, T_Storage>& data) {
  return o__gentype_eig(data);
} // eig

// user-callable: void eig ( std::complex<T> )
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void eig(const gMatrix<std::complex<T>, T_Storage1>& data,
         gMatrix<std::complex<T>, T_Storage2>& evectors,
         gMatrix<std::complex<T>, T_Storage3>& evalues) {
  o__gentype_eig(data, evectors, evalues);
} // eig

// Helper function that calls {c,z}geevx for the routines returning gEigenData
// eig of T returned as EigenData<gCppTraits<T>::T_Complex, gCppTraits<T>::T_Complex>
template<typename T, typename T_Storage>
gEigenData<typename gCppTraits<T>::T_Complex, typename gCppTraits<T>::T_Complex>
o__gentype_eig(const gMatrix<T, T_Storage>& data) {
  typedef typename gCppTraits<T>::T_Complex evect_type;
  typedef evect_type eval_type;
  // alloc matrices
  gMatrix<evect_type, ROW_DENSE> evectors;
  gMatrix<eval_type, ROW_DENSE> evalues;
  // call peer function
  o__gentype_eig(data, evectors, evalues);
  return gEigenData<evect_type, eval_type>(evectors, evalues);
} // o__gentype_eig

// Helper function that calls {c,z}geevx for the routines returning result in args.
// eig of T returned into matrices of gCppTraits<T>::T_Complex
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void o__gentype_eig(const gMatrix<T, T_Storage1>& data,
                    gMatrix<typename gCppTraits<T>::T_Complex, T_Storage2>& evectors,
                    gMatrix<typename gCppTraits<T>::T_Complex, T_Storage3>& evalues) {
#ifdef CHK_IDX
  if ( data.rows() != data.cols() ) {
    gErrorShape( "",
     "Eigen Analysis",
     "Matrix not square",
     data.rows(),
     data.cols() );
  } // if
#endif // CHK_IDX
  fortran::integer N = data.rows();
  typedef typename gCppTraits<T>::T_Real r_type;
  typedef typename gFortranTraits<r_type>::T_Complex T_Fortran;
  fortran::fortran_array<T_Fortran> a(N*N);
  // init fortran matrix with complex values. use helper. (this converts real input
  // to complex). This is done so that real AND complex input matrices may be handled
  // by cgeevx or zgeevx (complex Lapack routines). The real routines (sgeevx and dgeevx)
  // require some conversion and real/imag combining afterwards (that sux).
  o__copy_c_to_f(data, a);

  fortran::fortran_array<T_Fortran> vr(N*N);
  fortran::fortran_array<T_Fortran> w(N);

  // dispatch the Fortran arrays to the Lapack expert driver, {c,z}geevx (hooooo)
  o__geevx_dispatch(a, N, vr, w);
  // resize evector matrix and copy returned vctors
  evectors.resize(N, N);
  o__copy_f_to_c(vr, evectors);
  // resize evalue matrix, zero all elements and copy returned values onto diagonal
  evalues.resize(N, N);
  // JADP 17 Maart 2002 evalues = (T) 0;
  evalues.assignToAll( static_cast<T>(0) );
  for (T_Index rindx = 0; rindx < N; rindx++) {
    o__copy_f_to_c( w[rindx], evalues(rindx, rindx) );
  } // for
} // o__gentype_eig

/* #################### Singular Value Decomposition (SVD) ############  */

// user-callable void svd( float )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<float, T_Storage1>& A,
         gMatrix<float, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<float, T_Storage4>& VT) {
  o__gentype_svd(A, U, S, VT);
} // svd

// user-callable void svd( std::complex<float> )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<std::complex<float>, T_Storage1>& A,
         gMatrix<std::complex<float>, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<std::complex<float>, T_Storage4>& VT) {
  o__gentype_svd(A, U, S, VT);
} // svd

// user-callable void svd( double )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<double, T_Storage1>& A,
         gMatrix<double, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<double, T_Storage4>& VT) {
  o__gentype_svd(A, U, S, VT);
} // svd

// user-callable void svd( std::complex<double> )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<std::complex<double>, T_Storage1>& A,
         gMatrix<std::complex<double>, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<std::complex<double>, T_Storage4>& VT) {
  o__gentype_svd(A, U, S, VT);
} // svd

template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void o__gentype_svd(const gMatrix<T, T_Storage1>& A,
                    gMatrix<T, T_Storage2>& U,
                    gMatrix<typename gFortranTraits<T>::T_Real, T_Storage3>& S,
                    gMatrix<T, T_Storage4>& VT) {
  // some handy types
  typedef typename gFortranTraits<T>::T_Fortran T_Fortran;
  typedef typename gFortranTraits<T_Fortran>::T_Real T_Real;
  typedef typename gMatrix<T, ROW_DENSE>::T_Index T_Index;
  // declare locals
  const fortran::integer M = A.rows();
  const fortran::integer N = A.cols();
  const fortran::integer minMN = ( GLINEAR_o__MIN(M, N) );

  fortran::fortran_array<T_Fortran> a(M*N);
  fortran::fortran_array<T_Fortran> u(M*M);
  fortran::fortran_array<T_Fortran> vt(N*N);
  fortran::fortran_array<T_Real> s(minMN);
  // init Fortran matrix
  o__copy_c_to_f(A, a);
  // dispatch
  o__gesvd_dispatch(a, M, N, s, u, vt);
  // copy results back into C space
  U.resize(M, M);
  o__copy_f_to_c(u, U);
  VT.resize(N, N);
  o__copy_f_to_c(vt, VT);
  S.resize(minMN, minMN);
  // JADP 17Maart2002 S = static_cast<gFortranTraits<T>::T_Real>(0);
  S.assignToAll(static_cast<typename gFortranTraits<T>::T_Real>(0) );
  for (T_Index rindx = 0; rindx < minMN; rindx++) {
    o__copy_f_to_c( s[rindx], S(rindx, rindx) );
  } // for

} // o__gentype_svd

/* #################### Reduced Singular Value Decomposition (SVD) ############  */

// user-callable void reduced_svd( float )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<float, T_Storage1>& A,
         gMatrix<float, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<float, T_Storage4>& VT) {
  o__gentype_reduced_svd(A, U, S, VT);
} // reduced_svd

// user-callable void reduced_svd( std::complex<float> )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<std::complex<float>, T_Storage1>& A,
         gMatrix<std::complex<float>, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<std::complex<float>, T_Storage4>& VT) {
  o__gentype_reduced_svd(A, U, S, VT);
} // reduced_svd

// user-callable void reduced_svd( double )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<double, T_Storage1>& A,
         gMatrix<double, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<double, T_Storage4>& VT) {
  o__gentype_reduced_svd(A, U, S, VT);
} // reduced_svd

// user-callable void reduced_svd( std::complex<double> )
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<std::complex<double>, T_Storage1>& A,
         gMatrix<std::complex<double>, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<std::complex<double>, T_Storage4>& VT) {
  o__gentype_reduced_svd(A, U, S, VT);
} // reduced_svd

template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void o__gentype_reduced_svd(const gMatrix<T, T_Storage1>& A,
                    gMatrix<T, T_Storage2>& U,
                    gMatrix<typename gFortranTraits<T>::T_Real, T_Storage3>& S,
                    gMatrix<T, T_Storage4>& VT) {
  // some handy types
  typedef typename gFortranTraits<T>::T_Fortran T_Fortran;
  typedef typename gFortranTraits<T_Fortran>::T_Real T_Real;
  typedef typename gMatrix<T, ROW_DENSE>::T_Index T_Index;
  // declare locals
  const fortran::integer M = A.rows();
  const fortran::integer N = A.cols();
  const fortran::integer minMN = ( GLINEAR_o__MIN(M, N) );

  fortran::fortran_array<T_Fortran> a(M*N);
  fortran::fortran_array<T_Fortran> u(1);
  fortran::fortran_array<T_Fortran> vt(minMN*N);
  fortran::fortran_array<T_Real> s(minMN);
  // init Fortran matrix
  o__copy_c_to_f(A, a);
  // dispatch
  o__reduced_gesvd_dispatch(a, M, N, s, u, vt);
  // copy results back into C space
  VT.resize(minMN, N);
  o__copy_f_to_c(vt, VT);
  S.resize(minMN, minMN);
  for (T_Index rindx = 0; rindx < minMN; rindx++) {
    o__copy_f_to_c( s[rindx], S(rindx, rindx) );
  } // for

} // o__gentype_reduced_svd


GLINEAR_NAMESPACE_END
