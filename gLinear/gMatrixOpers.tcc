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
#include <algorithm>  // transform
#include <iostream>  // ostream, ios, endl, istream
#include <cstdlib>  // abs

GLINEAR_NAMESPACE_BEGIN

/*
typename gMatrix<T, T_Storage>::T_Dense
impulse(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
                 funcTransformVecInto< funcImpulse<T> >( funcImpulse<T>() ) );
  return result;
*/



// unary minus
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
operator-(const gMatrix<T, T_Storage>& arg) {
  typedef gMatrix<T, T_Storage> T_Matrix;
  typedef typename T_Matrix::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcMinus<T> >( funcMinus<T>() ) );
  return result;
}

// zeros
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> zeros(typename gMatrix<T, T_Storage>::T_Index rows,
          typename gMatrix<T, T_Storage>::T_Index cols) {
  gMatrix<T, T_Storage> tmp(rows, cols);
  // JADP 17 Maart 2002 tmp = 0;
  tmp.assignToAll(0);
  return tmp;
}

template<typename T>
gMatrix<T, ROW_DENSE> zeros(typename gMatrix<T, ROW_DENSE>::T_Index rows,
          typename gMatrix<T, ROW_DENSE>::T_Index cols) {
  return zeros<T, ROW_DENSE>(rows, cols);
}


// ones
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> ones(typename gMatrix<T, T_Storage>::T_Index rows,
         typename gMatrix<T, T_Storage>::T_Index cols) {
  gMatrix<T, T_Storage> tmp(rows, cols);
  // JADP 17 Maart 2002 tmp = 1;
  tmp.assignToAll(1);
  return tmp;
}

template<typename T>
gMatrix<T, ROW_DENSE> ones(typename gMatrix<T, ROW_DENSE>::T_Index rows,
         typename gMatrix<T, ROW_DENSE>::T_Index cols) {
  return ones<T, ROW_DENSE>(rows, cols);
}


// identity gMatrix
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> eye(typename gMatrix<T, T_Storage>::T_Index rows,
        typename gMatrix<T, T_Storage>::T_Index cols) {
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  if (cols == 0) {
    cols = rows;
  }
  gMatrix<T, T_Storage> tmp(rows, cols);
  // JADP 17 Maart 2002 tmp = 0;
  tmp.assignToAll(0);
  const T_Index mindim = ( (rows < cols) ? rows : cols );
  for (T_Index indx = 0; indx < mindim; indx++) {
    tmp(indx, indx) = 1;
  }
  return tmp;
}

template<typename T>
gMatrix<T, ROW_DENSE> eye(typename gMatrix<T, ROW_DENSE>::T_Index rows,
        typename gMatrix<T, ROW_DENSE>::T_Index cols) {
  return eye<T, ROW_DENSE>(rows, cols);
}

// diagonal matrix from row vector: interface -- forward to driver
template<typename T, typename T_VecStor>
gMatrix<T, ROW_DENSE> diag(const gRowVector<T, T_VecStor>& vec,
         int K) {
  return o__diagMatfromVec<gRowVector<T, T_VecStor>, ROW_DENSE>(vec, K);
}

// diagonal matrix from colvector: interface -- forward to driver
template<typename T, typename T_VecStor>
gMatrix<T, ROW_DENSE> diag(const gColVector<T, T_VecStor>& vec,
         int K) {
  return o__diagMatfromVec<gColVector<T, T_VecStor>, ROW_DENSE>(vec, K);
}

/*
template<typename T_Vector>
gMatrix<typename T_Vector::T_Value, ROW_DENSE> diag(const T_Vector& vec,
                                                    int K = 0) {
  return diag<T_Vector, ROW_DENSE>(vec, K);
}
*/


// diagonal matrix from vector: driver --- does the work
template<typename T_Vector, typename T_Storage>
gMatrix<typename T_Vector::T_Value, T_Storage>
o__diagMatfromVec(const T_Vector& vec,
      int K) {
  typedef typename T_Vector::T_Value T_Value;
  typedef gMatrix<T_Value, T_Storage> T_Matrix;
  typedef typename T_Matrix::T_Index T_Index;
  T_Index dim = vec.size();
  int aK = std::abs(K);
  T_Matrix tmp(dim + aK, dim + aK);
  // JADP 17 Maart 2002 tmp = 0;
  tmp.assignToAll(0);
  for (T_Index indx = 0; indx < dim; indx++) {
    if (K >= 0) {
      tmp(indx, indx + aK) = vec[indx];
    }
    else {
      tmp(indx + aK, indx) = vec[indx];
    }
  }
  return tmp;
}

// rowvector from matrix diagonal
template<typename T, typename T_Storage>
gRowVector<T, DENSE> diag(const gMatrix<T, T_Storage>& mat) {
  typedef gRowVector<T, DENSE> T_Vector;
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  T_Index dim = mat.rows();
  T_Vector tmp(dim);
  // JADP 17 Maart 2002 tmp = 0;
  tmp.assignToAll(0);
  for (T_Index indx = 0; indx < dim; indx++) {
    tmp[indx] = mat(indx, indx);
  }
  return tmp;
}

// general colvector - rowvector multiply: outerproduct
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
gMatrix<typename gNumPrec2<T1, T2>::T_Highest, ROW_DENSE>
operator*(const gColVector<T1, T_LHStor>& lhvec,
    const gRowVector<T2, T_RHStor>& rhvec) {
  typedef gMatrix<typename gNumPrec2<T1, T2>::T_Highest, ROW_DENSE> T_Matrix;
  typedef typename T_Matrix::T_Index T_Index;
  const T_Index rows = lhvec.size();
  const T_Index cols = rhvec.size();
  gMatrix<typename gNumPrec2<T1, T2>::T_Highest, ROW_DENSE> mat(rows, cols);
  typename gColVector<T1, T_LHStor>::T_ConstIterator cviter = lhvec.begin();
  const typename T_Matrix::T_ConstRowIterator mriter_end = mat.rowEnd();
  for (typename T_Matrix::T_RowIterator mriter = mat.rowBegin();
       mriter != mriter_end; ++mriter, ++cviter ) {
    *mriter = *cviter * rhvec;
  }
  return mat;
}

// general Matrix - Matrix equality operator
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
inline bool
operator==(const gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( ( lhmat.rows() != rhmat.rows() ) || ( lhmat.cols() != rhmat.cols() ) ) {
    gErrorSize( "",
    "operator==(const gMatrix<T1, T_LHStor>&, const gMatrix<T2, T_RHStor>&)",
    "Matrix sizes disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  const typename gMatrix<T1, T_LHStor>::T_ConstRowIterator lhiter_end = lhmat.rowEnd();
  typename gMatrix<T2, T_RHStor>::T_ConstRowIterator rhiter = rhmat.rowBegin();
  for (typename gMatrix<T1, T_LHStor>::T_ConstRowIterator lhiter = lhmat.rowBegin();
       lhiter != lhiter_end;
       lhiter++) {
    if ( ! (*lhiter == *rhiter++) ) {
      return false;
    } // if
  }
  return true;
}

// general Matrix - Matrix inequality operator
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
inline bool
operator!=(const gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat) {
  return !( operator==(lhmat, rhmat) );
}

// general gMatrix - gMatrix addition
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator+(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( ( lhmat.rows() != rhmat.rows() ) || ( lhmat.cols() != rhmat.cols() ) ) {
    gErrorSize( "",
    "operator+(const gMatrix<T1, T_LHStor>&, const gMatrix<T2, T_RHStor>&)",
    "Matrix sizes disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense T_Result;
  T_Result result( lhmat.rows(), lhmat.cols() );
  transformInto( lhmat.rowBegin(), lhmat.rowEnd(), rhmat.rowBegin(), result.rowBegin(),
     funcTransformVecInto< funcAdd<T1, T2> >( funcAdd<T1, T2>() ) );
  return result;
}


// general gMatrix - gMatrix subtraction
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator-(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( ( lhmat.rows() != rhmat.rows() ) || ( lhmat.cols() != rhmat.cols() ) ) {
    gErrorSize( "",
    "operator-(const gMatrix<T1, T_LHStor>&, const gMatrix<T2, T_RHStor>&)",
    "Matrix sizes disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense T_Result;
  T_Result result( lhmat.rows(), lhmat.cols() );
  transformInto( lhmat.rowBegin(), lhmat.rowEnd(), rhmat.rowBegin(), result.rowBegin(),
     funcTransformVecInto< funcSubtract<T1, T2> >( funcSubtract<T1, T2>() ) );
  return result;
}

// general gMatrix - gMatrix multiplication
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator*(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( lhmat.cols() != rhmat.rows() ) {
    gErrorSize( "",
    "operator*(const gMatrix<T, T_LHStor>&, const gMatrix<T, T_RHStor>&)",
    "Inner Matrix dimensions disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense T_Result;
  typedef typename T_Result::T_Index T_Index;
  T_Index rows, cols;
  rows = lhmat.rows();
  cols = rhmat.cols();
  T_Result result(rows, cols);

  for (T_Index r = 0; r < rows; r++) {
    for (T_Index c = 0; c < cols; c++) {
      result(r, c) = lhmat.row(r) * rhmat.col(c);
    }
  }
  return result;
}

// general gMatrix - gMatrix computed assignment: addition
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
inline gMatrix<T1, T_LHStor>&
operator+=(gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( ( lhmat.rows() != rhmat.rows() ) || ( lhmat.cols() != rhmat.cols() ) ) {
    gErrorSize( "",
    "operator+=(const gMatrix<T1, T_LHStor>&, const gMatrix<T2, T_RHStor>&)",
    "Matrix sizes disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  transformInto( lhmat.rowBegin(), lhmat.rowEnd(), rhmat.rowBegin(), lhmat.rowBegin(),
     funcTransformVecInto< funcAdd<T1, T2> >( funcAdd<T1, T2>() ) );
  return lhmat;
}

// general gMatrix - gMatrix computed assignment: subtraction
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
inline gMatrix<T1, T_LHStor>&
operator-=(gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat) {
#ifdef CHK_IDX
  if ( ( lhmat.rows() != rhmat.rows() ) || ( lhmat.cols() != rhmat.cols() ) ) {
    gErrorSize( "",
    "operator-=(const gMatrix<T, T_LHStor>&, const gMatrix<T, T_RHStor>&)",
    "Matrix sizes disagree",
    lhmat.rows(), lhmat.cols(), rhmat.rows(), rhmat.cols() );
  }
#endif // CHK_IDX
  transformInto( lhmat.rowBegin(), lhmat.rowEnd(), rhmat.rowBegin(), lhmat.rowBegin(),
     funcTransformVecInto< funcSubtract<T1, T2> >( funcSubtract<T1, T2>() ) );
  return lhmat;
}

// general gMatrix - colvector multiplication
template<typename T1, typename T2, typename T_MatStor, typename T_VecStor>
typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense
operator*(const gMatrix<T1, T_MatStor>& mat,
    const gColVector<T2, T_VecStor>& vec) {
#ifdef CHK_IDX
  if ( mat.cols() != vec.size() ) {
    gErrorSize( "",
    "operator*(const gMatrix<T1, T_MatStor>&, const gColVector<T2, T_VecStor>&)",
    "Matrix columns should equal vector size",
    mat.rows(), mat.cols(), vec.size() );
  }
#endif // CHK_IDX
  typedef typename gMatrix<T1, T_MatStor>::T_Row T_Row;
  typedef gColVector<T2, T_VecStor> T_Vector;
  typedef typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense T_Result;
  T_Result result( mat.rows() );
  std::transform( mat.rowBegin(), mat.rowEnd(), result.begin(), funcDot<T_Row, T_Vector>(vec) );
  return result;
}

// general rowvector - gMatrix multiplication
template<typename T1, typename T2, typename T_MatStor, typename T_VecStor>
typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense
// typename gRowVector<T, T_VecStor>::T_Dense
operator*(const gRowVector<T1, T_VecStor>& vec,
    const gMatrix<T2, T_MatStor>& mat) {
#ifdef CHK_IDX
  if ( mat.rows() != vec.size() ) {
    gErrorSize( "",
    "operator*(const gRowVector<T1, T_VecStor>&), const gMatrix<T2, T_MatStor>&)",
    "Matrix rows should equal vector size",
    mat.rows(), mat.cols(), vec.size() );
  }
#endif // CHK_IDX
  typedef typename gMatrix<T2, T_MatStor>::T_Col T_Col;
  typedef gRowVector<T1, T_VecStor> T_Vector;
  typedef typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense T_Result;
  T_Result result( mat.cols() );
  std::transform( mat.colBegin(), mat.colEnd(), result.begin(), funcDot<T_Col, T_Vector>(vec) );
  return result;
}

// general scalar - gMatrix addition
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator+(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcAdd<T2, T1> >( funcAdd<T2, T1>(scalar) ) );
  return result;
}

// general scalar - gMatrix subtraction
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator-(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcSubtract<T2, T1, 1> >( funcSubtract<T2, T1, 1>(scalar) ) );
  return result;
}

// general scalar - gMatrix multiplication
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator*(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcMultiply<T2, T1> >( funcMultiply<T2, T1>(scalar) ) );
  return result;
}

// general scalar - gMatrix division
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator/(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcDivide<T2, T1, 1> >( funcDivide<T2, T1, 1>(scalar) ) );
  return result;
}

// general gMatrix - scalar addition
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator+(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar) {
  return scalar + mat;
}

// general gMatrix - scalar subtraction
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator-(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcSubtract<T1, T2> >( funcSubtract<T1, T2>(scalar) ) );
  return result;
}

// general gMatrix - scalar multiplication
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator*(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar) {
  return scalar * mat;
}

// general gMatrix - scalar division
template<typename T1, typename T2, typename T_Storage>
inline typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator/(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar) {
  typedef typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense T_Result;
  T_Result result( mat.rows(), mat.cols() );
  transformInto( mat.rowBegin(), mat.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcDivide<T1, T2> >( funcDivide<T1, T2>(scalar) ) );
  return result;
}

// general gMatrix - scalar computed assignment: addition
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>&
operator+=(gMatrix<T1, T_Storage>& mat,
     const T2& scalar) {
  transformInto( mat.rowBegin(), mat.rowEnd(), mat.rowBegin(),
     funcTransformVecInto< funcAdd<T1, T2> >( funcAdd<T1, T2>(scalar) ) );
  return mat;
}

// general gMatrix - scalar computed assignment: subtraction
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>&
operator-=(gMatrix<T1, T_Storage>& mat,
     const T2& scalar) {
  transformInto( mat.rowBegin(), mat.rowEnd(), mat.rowBegin(),
     funcTransformVecInto< funcSubtract<T1, T2> >( funcSubtract<T1, T2>(scalar) ) );
  return mat;
}

// general gMatrix - scalar computed assignment: multiplication
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>&
operator*=(gMatrix<T1, T_Storage>& mat,
     const T2& scalar) {
  transformInto( mat.rowBegin(), mat.rowEnd(), mat.rowBegin(),
     funcTransformVecInto< funcMultiply<T1, T2> >( funcMultiply<T1, T2>(scalar) ) );
  return mat;
}

// general gMatrix - scalar computed assignment: division
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>&
operator/=(gMatrix<T1, T_Storage>& mat,
     const T2& scalar) {
  transformInto( mat.rowBegin(), mat.rowEnd(), mat.rowBegin(),
     funcTransformVecInto< funcDivide<T1, T2> >( funcDivide<T1, T2>(scalar) ) );
  return mat;
}


#if 0
#define GLINEAR_o__apply_gMatrix_transform(src, dst, func, r, c)  \
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;  \
  for (T_Index rindx = 0; rindx < r; rindx++) {      \
    for (T_Index cindx = 0; cindx < c; cindx++) {    \
      dst(rindx, cindx) = func( src(rindx,cindx) );    \
    }                \
  }
#endif // 0

// return gMatrix containing complex conjugate of elements of arg.
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
conj(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcConj<T> >( funcConj<T>() ) );
  return result;
}

// return gMatrix containing real components of elements of arg.
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
    typename gMatrixStorTraits<T_Storage>::T_Dense>
real(const gMatrix<T, T_Storage>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gMatrix<T_Real, typename gMatrixStorTraits<T_Storage>::T_Dense> T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcReal<T> >( funcReal<T>() ) );
  return result;
}

// return gMatrix containing imaginary components of elements of arg.
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
  typename gMatrixStorTraits<T_Storage>::T_Dense>
imag(const gMatrix<T, T_Storage>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gMatrix<T_Real, typename gMatrixStorTraits<T_Storage>::T_Dense> T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcImag<T> >( funcImag<T>() ) );
  return result;
}

// absolute value
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
  typename gMatrixStorTraits<T_Storage>::T_Dense>
abs(const gMatrix<T, T_Storage>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gMatrix<T_Real, typename gMatrixStorTraits<T_Storage>::T_Dense> T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcAbs<T> >( funcAbs<T>() ) );
  return result;
}

// argument (angle)
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
  typename gMatrixStorTraits<T_Storage>::T_Dense>
arg(const gMatrix<T, T_Storage>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gMatrix<T_Real, typename gMatrixStorTraits<T_Storage>::T_Dense> T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcArg<T> >( funcArg<T>() ) );
  return result;
}

// norm (abs ^ 2)
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
  typename gMatrixStorTraits<T_Storage>::T_Dense>
norm(const gMatrix<T, T_Storage>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gMatrix<T_Real, typename gMatrixStorTraits<T_Storage>::T_Dense> T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcNorm<T> >( funcNorm<T>() ) );
  return result;
}

// cosine
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
cos(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcCos<T> >( funcCos<T>() ) );
  return result;
}

// hyperbolic cosine
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
cosh(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcCosh<T> >( funcCosh<T>() ) );
  return result;
}

// sine
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sin(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcSin<T> >( funcSin<T>() ) );
  return result;
}

// sinh
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sinh(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcSinh<T> >( funcSinh<T>() ) );
  return result;
}

// exp
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
exp(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcExp<T> >( funcExp<T>() ) );
  return result;
}

// log
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
log(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcLog<T> >( funcLog<T>() ) );
  return result;
}

// sqrt
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sqrt(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcSqrt<T> >( funcSqrt<T>() ) );
  return result;
}

// step function
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
step(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcStep<T> >( funcStep<T>() ) );
  return result;
}

// impulse function
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
impulse(const gMatrix<T, T_Storage>& arg) {
  typedef typename gMatrix<T, T_Storage>::T_Dense T_Result;
  T_Result result( arg.rows(), arg.cols() );
  transformInto( arg.rowBegin(), arg.rowEnd(), result.rowBegin(),
     funcTransformVecInto< funcImpulse<T> >( funcImpulse<T>() ) );
  return result;
}

// write dims and contents to stream - human readable
template<typename T, typename T_Storage>
std::ostream& operator<<(std::ostream& os,
       const gMatrix<T, T_Storage>& mat) {
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  const T_Index nRows = mat.rows();
  const T_Index nCols = mat.cols();
  // long oldFlags; JADP 3May2003
  std::ios::fmtflags oldFlags;

  os << nRows << " " << nCols << std::endl;
  oldFlags = os.flags(std::ios::scientific);
  for (int rindx = 0; rindx < nRows; rindx++) {
    for (int cindx = 0; cindx < nCols; cindx++) {
      os << mat(rindx, cindx) << " ";
    }
    os << std::endl;
  }
  os.flags(oldFlags);
  return os;
}

// read dims and contents from stream - human readable
template<typename T, typename T_Storage>
std::istream& operator>>(std::istream& is,
       gMatrix<T, T_Storage>& mat) {
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  T_Index nRows;
  T_Index nCols;
  is >> nRows >> nCols;
  mat.resize(nRows, nCols);
  for (int rindx = 0; rindx < nRows; rindx++) {
    for (int cindx = 0; cindx < nCols; cindx++) {
      is >> mat(rindx, cindx);
    }
  }
  return is;
}

// write dims and contents to stream - binary
// NOTE: This is platform DEPENDENT, as T_Index is normally typedef'ed to
// size_t which is big enough to hold the size of a pointer (which differs
// between 32- and 64-bit archis). I warned you.
template<typename T, typename T_Storage>
std::ostream& global_write(std::ostream& os,
               const gMatrix<T, T_Storage>& mat) {
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  const T_Index nRows = mat.rows();
  const T_Index nCols = mat.cols();

  os.write( reinterpret_cast<const char*>(&nRows), sizeof(nRows) );
  os.write( reinterpret_cast<const char*>(&nCols), sizeof(nCols) );
  for (int rindx = 0; rindx < nRows; rindx++) {
    for (int cindx = 0; cindx < nCols; cindx++) {
      os.write( reinterpret_cast<const char*>( &mat(rindx, cindx) ), sizeof(T) );
    }
  }
  return os;
}

// read dims and contents from stream - binary
// NOTE: This is platform DEPENDENT, as T_Index is normally typedef'ed to
// size_t which is big enough to hold the size of a pointer (which differs
// between 32- and 64-bit archis). I warned you.
template<class T, typename T_Storage>
std::istream& global_read(std::istream& is,
              gMatrix<T, T_Storage>& mat) {
  typedef typename gMatrix<T, T_Storage>::T_Index T_Index;
  T_Index nRows;
  T_Index nCols;

  is.read( reinterpret_cast<char *>(&nRows), sizeof(nRows) );
  is.read( reinterpret_cast<char *>(&nCols), sizeof(nCols) );
  mat.resize(nRows, nCols);
  for (int rindx = 0; rindx < nRows; rindx++) {
    for (int cindx = 0; cindx < nCols; cindx++) {
      is.read( reinterpret_cast<char*>( &mat(rindx, cindx) ), sizeof(T) );
    }
  }
  return is;
}

GLINEAR_NAMESPACE_END
