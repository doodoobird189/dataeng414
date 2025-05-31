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

#ifndef GMATRIX_OPERS_H
#define GMATRIX_OPERS_H

// glinear headers
#include "gLinearMacros.h"
#include "gRowVector.h"
#include "gColVector.h"
#include "gMatrix.h"
#include "gLinearAlgo.h"
#include "gLinearFunctors.h"
#include "gCppTraits.h"
#include "gMatrixTraits.h"

// standard headers
#include <iosfwd>  // ostream, istream

GLINEAR_NAMESPACE_BEGIN


/**@name Unary Matrix operators
 */
//@{
/**
 * gMatrix<T, T_Stor> unary minus.
 * Return matrix with operator-(T) applied to each element.
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense operator-(const gMatrix<T, T_Storage>& arg);
//@}

/**@name Some Handy Matrices
 */
//@{
/**
 * Matrix filled with 0's. \\
 * Consruct a [rows x cols] matrix filled with zeros.
 * @param rows Number of rows
 * @param cols Number of columns
 * @return [rows x cols] matrix filled with zeros.
 */
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> zeros(typename gMatrix<T, T_Storage>::T_Index rows,
          typename gMatrix<T, T_Storage>::T_Index cols);
template<typename T>
gMatrix<T, ROW_DENSE> zeros(typename gMatrix<T, ROW_DENSE>::T_Index rows,
          typename gMatrix<T, ROW_DENSE>::T_Index cols);

/**
 * Matrix filled with 1's. \\
 * Consruct a [rows x cols] matrix filled with ones.
 * @param rows Number of rows
 * @param cols Number of columns
 * @return [rows x cols] matrix filled with ones.
 */
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> ones(typename gMatrix<T, T_Storage>::T_Index rows,
                           typename gMatrix<T, T_Storage>::T_Index cols);
template<typename T>
gMatrix<T, ROW_DENSE> ones(typename gMatrix<T, ROW_DENSE>::T_Index rows,
                           typename gMatrix<T, ROW_DENSE>::T_Index cols);

/**
 * Identity matrix. \\
 * Construct a [rows x cols] identity matrix.
 * @param rows Number of rows
 * @param cols Number of columns (optional). if (cols == 1) cols = rows;
 * @return [rows x cols] identity matrix.
 */
template<typename T, typename T_Storage>
gMatrix<T, T_Storage> eye(typename gMatrix<T, T_Storage>::T_Index rows,
                          typename gMatrix<T, T_Storage>::T_Index cols = 0);
template<typename T>
gMatrix<T, ROW_DENSE> eye(typename gMatrix<T, ROW_DENSE>::T_Index rows,
                          typename gMatrix<T, ROW_DENSE>::T_Index cols = 0);

//@}

/**@name Vector - Matrix conversions
 */
//@{

// driver function
template<typename T_Vector, typename T_Storage>
gMatrix<typename T_Vector::T_Value, T_Storage> o__diagMatfromVec(const T_Vector& vec, int K);

/**
 * Diagonal from RowVector. \\
 * Construct an [N x N] matrix from an RowVector.
 * The vector elements are the elements on the K'th matrix diagonal
 * and all the other elements are zero. The main diagonal is the 0'th
 * diagonal, the one below the main is -1 and the one above is 1.
 */
template<typename T, typename T_VecStor>
gMatrix<T, ROW_DENSE> diag(const gRowVector<T, T_VecStor>& vec,
                           int K = 0);

/**
 * Diagonal from ColVector. \\
 * Construct an [N x N] matrix from an RowVector.
 * The vector elements are the elements on the K'th matrix diagonal
 * and all the other elements are zero. The main diagonal is the 0'th
 * diagonal, the one below the main is -1 and the one above is 1.
 */
template<typename T, typename T_VecStor>
gMatrix<T, ROW_DENSE> diag(const gColVector<T, T_VecStor>& vec,
                           int K = 0);

/*
template<typename T_Vector>
gMatrix<typename T_Vector::T_Value, ROW_DENSE> diag(const T_Vector& vec,
                int K = 0);
*/

/**
 * Rowvector from diagonal. \\
 * Construct an N-length rowvector from the matrix diagonal.
 */
template<typename T, typename T_Storage>
gRowVector<T, DENSE> diag(const gMatrix<T, T_Storage>& mat);

// JADP 3May2003 the above was as below, but out of sync with the tcc implementation.
// Changed to DENSE and removed the K parameter believed to be an error
// template<typename T, typename T_Storage>
// gRowVector<T, ROW_DENSE> diag(const gMatrix<T, T_Storage>& mat,
//         int K = 0);

/**@name Vector - Vector
 */
//@{
/**
 * Vector - Vector outer product
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
gMatrix<typename gNumPrec2<T1, T2>::T_Highest, ROW_DENSE>
operator*(const gColVector<T1, T_LHStor>& lhvec,
    const gRowVector<T2, T_RHStor>& rhvec);
//@}

/**@name Matrix - Matrix
 */
//@{
/**
 * General Matrix - Matrix equality operator
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
bool operator==(const gMatrix<T1, T_LHStor>& lhmat,
                const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General Matrix - Matrix inequality operator
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
bool operator!=(const gMatrix<T1, T_LHStor>& lhmat,
                const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General Matrix - Matrix addition operator
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator+(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General Matrix - Matrix subtraction operator
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator-(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General Matrix - Matrix multiplication operator
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_LHStor>::T_Dense
operator*(const gMatrix<T1, T_LHStor>& lhmat,
    const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General computed assignment: addition
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
gMatrix<T1, T_LHStor>&
operator+=(gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat);
/**
 * General computed assignment: subtraction
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
gMatrix<T1, T_LHStor>&
operator-=(gMatrix<T1, T_LHStor>& lhmat,
     const gMatrix<T2, T_RHStor>& rhmat);
//@}

/**@name Matrix - Vector
 */
//@{
/**
 * General Matrix - ColVector multiplication
 */
template<typename T1, typename T2, typename T_MatStor, typename T_VecStor>
typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense
operator*(const gMatrix<T1, T_MatStor>& mat,
    const gColVector<T2, T_VecStor>& vec);
/**
 * General RowVector - Matrix multiplication
 */
template<typename T1, typename T2, typename T_MatStor, typename T_VecStor>
typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, T_VecStor>::T_Dense
operator*(const gRowVector<T1, T_VecStor>& vec,
    const gMatrix<T2, T_MatStor>& mat);
//@}

/**@name Scalar - Matrix
*/
//@{
/**
 * General Scalar - Matrix addition
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator+(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat);

/**
 * General Scalar - Matrix subtraction
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator-(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat);
/**
 * General Scalar - Matrix multiplication
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator*(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat);

/**
 * General Scalar - Matrix division
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator/(const T1& scalar,
    const gMatrix<T2, T_Storage>& mat);
//@}

/**@name Matrix - Scalar
 */
//@{
/**
 * General Matrix - Scalar addition
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator+(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar);
/**
 * General Matrix - Scalar subtraction
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator-(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar);
/**
 * General Matrix - Scalar multiplication
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator*(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar);
/**
 * General Matrix - Scalar division
 */
template<typename T1, typename T2, typename T_Storage>
typename gMatrix<typename gNumPrec2<T1, T2>::T_Highest, T_Storage>::T_Dense
operator/(const gMatrix<T1, T_Storage>& mat,
    const T2& scalar);

/**
 * General Matrix - Scalar computed assignment: addition
 */
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>&
operator+=(gMatrix<T1, T_Storage>& mat,
     const T2& scalar);

/**
 * General Matrix - Scalar computed assignment: subtraction
 */
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>& operator-=(gMatrix<T1, T_Storage>& mat,
                                   const T2& scalar);
/**
 * General Matrix - Scalar computed assignment: multiplication
 */
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>& operator*=(gMatrix<T1, T_Storage>& mat,
                                   const T2& scalar);

/**
 * General Matrix - Scalar computed assignment: division
 */
template<typename T1, typename T2, typename T_Storage>
gMatrix<T1, T_Storage>& operator/=(gMatrix<T1, T_Storage>& mat,
                                   const T2& scalar);

//@}


/*
template<typename T1, typename T_Storage1, typename T2, typename T_Storage2>
void o__apply_matrix_transform( const gMatrix<T1, stortype1>& mat1,
                                gMatrix<T2, stortype2>& mat2,
                                T2(*transform)(const T1&) );

*/

/**@name Functions
 */
//@{

/**
 * Return matrix containing the complex conugate of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense conj(const gMatrix<T, T_Storage>& arg);

/**
 * Return matrix containing the real components of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
        typename gMatrixStorTraits<T_Storage>::T_Dense>
real(const gMatrix<T, T_Storage>& arg);

/**
 * Return matrix containing the imaginary components of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
    typename gMatrixStorTraits<T_Storage>::T_Dense>
imag(const gMatrix<T, T_Storage>& arg);

/**
 * Return matrix containing the absolute value of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
        typename gMatrixStorTraits<T_Storage>::T_Dense>
abs(const gMatrix<T, T_Storage>& arg);


/**
 * Return matrix containing the argument of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
    typename gMatrixStorTraits<T_Storage>::T_Dense>
arg(const gMatrix<T, T_Storage>& arg);

/**
 * Return matrix containing the norm of the
 * argument matrix elements
 */
template<typename T, typename T_Storage>
gMatrix<typename gCppTraits<T>::T_Real,
    typename gMatrixStorTraits<T_Storage>::T_Dense>
norm(const gMatrix<T, T_Storage>& arg);

/**
 * Cosine
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
cos(const gMatrix<T, T_Storage>& arg);

/**
 * Hyperbolic Cosine
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
cosh(const gMatrix<T, T_Storage>& arg);

/**
 * Sine
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sin(const gMatrix<T, T_Storage>& arg);

/**
 * Hyperbolic Sine
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sinh(const gMatrix<T, T_Storage>& arg);

/**
 * Exponent
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
exp(const gMatrix<T, T_Storage>& arg);

/**
 * Natural logarithm
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
log(const gMatrix<T, T_Storage>& arg);

/**
 * Square root
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
sqrt(const gMatrix<T, T_Storage>& arg);

/**
 * Unit step function
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
step(const gMatrix<T, T_Storage>& mat);

/**
 * Impulse (dirac-delta) function
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense
impulse(const gMatrix<T, T_Storage>& mat);

//@}

/**@name Human-readable Stream I/O
 */
//@{
/**
 * Write matrix length, followed by contents, to ostream
 */
template<typename T, typename T_Storage>
std::ostream& operator<<(std::ostream& os, const gMatrix<T, T_Storage>& mat);

/**
 * Read matrix length, followed by contents, from istream
 */
template<class T, typename T_Storage>
std::istream& operator>>(std::istream& is, gMatrix<T, T_Storage>& mat);

//@}

/**@name Binary Stream I/O
 */
//@{
/**
 * Write matrix length, followed by contents, to ostream (raw binary)
 */
template<typename T, typename T_Storage>
std::ostream& global_write(std::ostream& os, const gMatrix<T, T_Storage>& mat);

/**
 * Read matrix length, followed by contents, from istream (raw binary)
 */
template<typename T, typename T_Storage>
std::istream& global_read(std::istream& is, gMatrix<T, T_Storage>& mat);

//@}

GLINEAR_NAMESPACE_END

#include "gMatrixOpers.tcc"

#endif // GMATRIX_OPERS_H
