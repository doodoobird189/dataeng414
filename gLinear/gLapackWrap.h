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

#ifndef GLAPACK_WRAP_H
#define GLAPACK_WRAP_H

// glinear headers
#include "gMatrix.h"
#include "gFortranTraits.h"
#include "gCppTraits.h"
#include "gLapackHeaders.h"
#include "gRowVector.h"

// standard headers
#include <complex>  // complex

GLINEAR_NAMESPACE_BEGIN

/**@name Language Support
 */
//@{

#if 0
/**
 * Class doing some handy typename lookups relating the type of the
 * data in the matrix to the type of the eigenvectors and eigenvalues. \\
 * Each EigenTraits<T> instance consists contains the following typedefs: \\\\
 * \begin{description}
 * \item[T] The type T
 * \item[real_type] The related real-valued numerical type
 * \item[complex_type] The related complex-valued numerical type
 * \item[evalue_type] Type of the elements of the Eigenvalue matrix
 * \item[evector_type] Type of the elements of the Eigenvector matrix
 * \end{description}
 */
template<typename T>
struct EigenTraits {
  typedef typename T::data_type data_type;
  typedef typename T::real_type real_type;
  typedef typename T::complex_type complex_type;
  typedef typename T::evalue_type evalue_type;
  typedef typename T::evector_type evector_type;
}; // EigenTraits

/**
 * EigenTraits for Matrix of {\em float}.
 */
template<>
struct EigenTraits<float> {
  /**
   * Data type
   */
  typedef float data_type;
  /**
   * Related real-valued type
   */
  typedef float real_type;
  /**
   * Related complex-valued type
   */
  typedef std::complex<float> complex_type;
  /**
   * Related Eigenvalue element type
   */
  typedef complex_type evalue_type;
  /**
   * Related Eigenvector element
   */
  typedef complex_type evector_type;
}; // EigenTraits

/**
 * EigenTraits for Matrix of {\em complex<float> }
 */
template<>
struct EigenTraits< std::complex<float> > {
  /**
   * Data type
   */
  typedef std::complex<float> data_type;
  /**
   * Related real-valued type
   */
  typedef float real_type;
  /**
   * Related complex-valued type
   */
  typedef std::complex<float> complex_type;
  /**
   * Related Eigenvalue element type
   */
  typedef complex_type evalue_type;
  /**
   * Related Eigenvector element type
   */
  typedef complex_type evector_type;
}; // EigenTraits

/**
 * EigenTraits for Matrix of {\em double}
 */
template<>
struct EigenTraits<double> {
  /**
   * Data type
   */
  typedef double data_type;
  /**
   * Related real-valued type
   */
  typedef double real_type;
  /**
   * Related complex-valued type
   */
  typedef std::complex<double> complex_type;
  /**
   * Related Eigenvalue element type
   */
  typedef complex_type evalue_type;
  /**
   * Related Eigenvector element type
   */
  typedef complex_type evector_type;
}; // EigenTraits

/**
 * EigenTraits for Matrix of {\em complex<double>}
 */
template<>
struct EigenTraits< std::complex<double> > {
  /**
   * Data type
   */
  typedef std::complex<double> data_type;
  /**
   * Related real-valued type
   */
  typedef double real_type;
  /**
   * Related complex-valued type
   */
  typedef std::complex<double> complex_type;
  /**
   * Related Eigenvalue element type
   */
  typedef complex_type evalue_type;
  /**
   * Related Eigenvector element type
   */
  typedef complex_type evector_type;
}; // EigenTraits
#endif // 0

/**
 * Class Containing an Eigenvalue and Eigenvector Matrix pair
 */
template<typename TVect, typename TVal = TVect>
class gEigenData {
public:
  /**
   * Type of the elements in the Eigenvector Matrix
   */
  typedef TVect T_EVector;
  /**
   * Type of the elements in the Eigenvalue Matrix.
   */
  typedef TVal T_EValue;
private:
  /**
   * [N x N] Eigenvector matrix (ROW_DENSE)
   */
  gMatrix<T_EVector, ROW_DENSE> eVectors;
  /**
   * [N x N] Eigenvalue matrix (ROW_DENSE)
   */
  gMatrix<T_EValue, ROW_DENSE> eValues;
public:
  /**
   * Construct and init with supplied matrices.
   */
  gEigenData(const gMatrix<T_EVector, ROW_DENSE>& newvect,
             const gMatrix<T_EValue, ROW_DENSE>& newval);

  /**
   * Assign to other EigenData object.
   * The assignment operator is invoked for the eigenvector and -value
   * matrices, thereby performing a deep copy.
   */
  template<typename T_srcVec, typename T_srcVal>
  gEigenData& operator=(const gEigenData<T_srcVec, T_srcVal>& src);

  /**
   * Return refererence to the Eigenvector matrix.
   * @return Reference to the Eigenvector matrix.
   */
  gMatrix<T_EVector, ROW_DENSE>& vectors(void);

  /**
   * Return const refererence to the Eigenvector matrix.
   * @return Const reference to the Eigenvector matrix.
   */
  const gMatrix<T_EVector, ROW_DENSE>& vectors(void) const;

  /**
   * Return refererence to the Eigenvalue matrix.
   * @return Reference to the Eigenvalue matrix.
   */
  gMatrix<T_EValue, ROW_DENSE>& values(void);

  /**
   * Return const refererence to the Eigenvalue matrix.
   * @return Const reference to the Eigenvalue matrix.
   */
  const gMatrix<T_EValue, ROW_DENSE>& values(void) const;

}; // gEigenData

//@}


/**@name Matrix Norms and Condition Numbers
 */
//@{

/**
 * Compute and return the P-norm of matrix A.
 * @param A Argument matrix
 * @param P Norm to calculate (default = 2).
 */
template<typename T, typename T_Storage>
T matrixNorm(const gMatrix<T, T_Storage>& A, int P = 2);

/**
 * Compute and return the Matrix 2-norm (largest singular value) of A. \\
 *  norm2(A) = max(svd(A)).
 * @param A Argument matrix
 * @return Matrix 2-norm of A
 */
template<typename T, typename T_Storage>
T matrixNorm2(const gMatrix<T, T_Storage>& A);

/**
 * Compute and return the P-norm condition number of A.
 * @param A Argument matrix
 * @param P Norm
 */
template<typename T, typename T_Storage>
T cond(const gMatrix<T, T_Storage>& A, int P = 2);

/**
 * Compute and return the 2-norm condition number of A. \\
 * cond2(A) = max(svd(A)) / min(svd(A))
 * @param A Argument matrix
 * @param P Norm
 */
template<typename T, typename T_Storage>
T cond2(const gMatrix<T, T_Storage>& A);

//@}


/**@name Matrix Inversion
 */
//@{

/**
 * Return inverse of Matrix as another Matrix. \\
 * The inverse of the argument Matrix is returned
 * as a dense Matrix of the same orientation as the
 * argument matrix.
 * @param A input matrix
 * @return Inverse of A
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense inv(const gMatrix<T, T_Storage>& A,
                                            int& fail);
/**
 * Return inverse of Matrix as another Matrix and determinant
 * in supplied argument. \\
 * The inverse of the argument Matrix is returned
 * as a dense Matrix of the same orientation as the
 * argument matrix. Another supplied argument is returned equal
 * to the determinant of the inverse of the argument matrix.
 * @param A input matrix
 * @param det Contains determinant of the {\bf inverse} on function return.
 * @return Inverse of A
 */
template<typename T, typename T_Storage>
typename gMatrix<T, T_Storage>::T_Dense inv(const gMatrix<T, T_Storage>& A,
                                            T& det,
                                            int& fail);

/**
 * Return inverse of Matrix in an argument Matrix. \\
 * The inverse of the argument Matrix is returned
 * into the second Matrix argument
 * @param A input matrix
 * @param invA Contains the inverse of A on return.
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void inv(const gMatrix<T, T_Storage1>& A,
         gMatrix<T, T_Storage2>& invA,
         int& fail);
/**
 * Return inverse of Matrix and determinant of inverse
 * in supplied arguments. \\
 * @param A input matrix
 * @param invA Contains the inverse of A on return.
 * @param det Contains determinant of {\bf invA} on function return.
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void inv(const gMatrix<T, T_Storage1>& A,
         gMatrix<T, T_Storage2>& invA,
         T& det,
         int& fail);
/**
 * Return the determinant of the matrix
 * @param A input matrix
 * @return Determinant of A
 */
template<typename T, typename T_Storage>
  T det(const gMatrix<T, T_Storage>& A,
        int& fail);

/**
 * Returns the cholesky decomposition of argument matrix
 * A is real and symmetric, the lower triangle is used.
 * Result is returned in a lower triangle matrix L. A = L*L'
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void cholesky_factor(const gMatrix<T, T_Storage1>& A,
                     gMatrix<T, T_Storage2>& L,
                     int& fail);

/**
 * Returns the cholesky decomposition of argument matrix
 * A is real and symmetric, the lower triangle is used.
 * Result is returned in a lower triangle matrix L. A = L*L'
 *
 * The result is kept in the provided Fortran structure to be
 * re-used by the Cholesky solvers that accept a precomputed
 * Cholesky factor.
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1>
void cholesky_factor(const gMatrix<T, T_Storage1>& A,
                     fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                     int& fail);

/**
 * Return the LU decomposition and permutation vector
 * of argument matrix
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void lu_factor(const gMatrix<T, T_Storage1>& A,
               gMatrix<T, T_Storage2>& LU,
               gRowVector<int>& P);

/**
 * apply permutation to matrix
 */
template<typename T, typename T_Storage1>
inline void lu_permute(gMatrix<T, T_Storage1>& A,
                       const gRowVector<int>& P);

/**
 * Linear Systems:
 * Param 'trans' specifies the form of the system of equations:
 *  'N':  A * X = B  (No transpose)
 *  'T':  A'* X = B  (Transpose)
 *  'C':  A'* X = B  (Conjugate transpose = Transpose)
 */

/**
 * Solve a Linear System Ax = b.
 * A is real and symmetric, the lower triangle is used.
 * Rhs (b) and solution (x) are column vectors
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void cholesky_solve(const gMatrix<T, T_Storage1>& A,
                    const gColVector<T, T_Storage2>& b,
                    gColVector<T, T_Storage3>& x,
                    int& fail);


/**
 * Solve a Linear System Ax = b with a precomputed Cholesky factorization L of A.
 * L is real and symmetric, the lower triangle is used.
 * Rhs (b) and solution (x) are column vectors
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1, typename T_Storage2>
void cholesky_solve(fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                    const gColVector<T, T_Storage1>& b,
                    gColVector<T, T_Storage2>& x,
                    int& fail);

/**
 * Solve multiple Linear Systems AX = B.
 * A is real and symmetric, the lower triangle is used.
 * Each column in B represents the rhs of a linear system with the
 * corresponding column in X its solution.
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void cholesky_solve(const gMatrix<T, T_Storage1>& A,
                    const gMatrix<T, T_Storage2>& B,
                    gMatrix<T, T_Storage3>& X,
                    int& fail);

/**
 * Solve multiple Linear Systems AX = B with a precomputed Cholesky factorization L of A.
 * L is real and symmetric, the lower triangle is used.
 * Each column in B represents the rhs of a linear system with the
 * corresponding column in X its solution.
 * fail is 0 if successful
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void cholesky_solve(fortran::fortran_array<typename gFortranTraits<T>::T_Fortran>& L,
                    const gMatrix<T, T_Storage2>& B,
                    gMatrix<T, T_Storage3>& X,
                    int& fail);

/**
 * Solve a Linear System.  Rhs (b) and solution (x) are column vectors
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& A,
              const gColVector<T, T_Storage2>& b,
              gColVector<T, T_Storage3>& x,
              int& fail,
              const char trans = 'N');

/**
 * Solve multiple Linear Systems.
 * Each column in b represents the rhs of a linear system with the
 * corresponding column in x its solution.
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& A,
              const gMatrix<T, T_Storage2>& b,
              gMatrix<T, T_Storage3>& x,
              int& fail,
              const char trans = 'N');

/**
 * Solve a Linear System with provided LU factorization.
 * Rhs (b) and solution (x) are column vectors
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& LU,
              const gRowVector<int>& P,
              const gColVector<T, T_Storage2>& b,
              gColVector<T, T_Storage3>& x,
              int& fail,
              const char trans = 'N');


/**
 * Solve multiple Linear Systems with provided LU factorization
 * Each column in b represents the rhs of a linear system with the
 * corresponding column in x its solution.
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void lu_solve(const gMatrix<T, T_Storage1>& LU,
              const gRowVector<int>& P,
              const gMatrix<T, T_Storage2>& b,
              gMatrix<T, T_Storage3>& x,
              int& fail,
              const char trans = 'N');


//@}

/**@name Eigenvalue decomposition
 */
//@{
/**
 * Return Eigenvectors and -values of a Matrix of {\bf T} as an
 * {\bf EigenData< complex<T>, complex<T> >} object.
 * @param data [N x N] Matrix of {\em T} containing data vectors.
 * @return EigenData object containing:
 * \begin{itemize}
 * \item [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvectors of {\em data} as columns.
 * \item [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvalues of {\em data} on the diagonal.
 * \end{itemize}
 */
template<typename T, typename T_Storage>
gEigenData< std::complex<T>, std::complex<T> > eig(const gMatrix<T, T_Storage>& data);
/**
 * Return Eigenvectors and -values of a Matrix of {\bf T} into the
 * supplied argument matrices of {\bf complex<T> }.
 * @param data [N x N] Matrix of {\em T} containing data vectors.
 * @param evalues [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvectors as columns on exit. (Ignored on entry).
 * @param evectors [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvalues on the diagonal on exit. (Ignored on entry).
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void eig(const gMatrix<T, T_Storage1>& data,
         gMatrix<std::complex<T>, T_Storage2>& evectors,
         gMatrix<std::complex<T>, T_Storage3>& evalues);
/**
 * Return Eigenvectors and -values of a Matrix of {\bf complex<T>} as an
 * {\bf EigenData< complex<T>, complex<T> >} object.
 * @param data [N x N] Matrix of {\em complex<T>} containing data vectors.
 * @return EigenData object containing:
 * \begin{itemize}
 * \item [N x N] Matrix of {\em complex<T>} containing the
 *  eigenvectors of {\em data} as columns.
 * \item [N x N] Matrix of {\em complex<T>} containing the
 *  eigenvalues of {\em data} on the diagonal.
 * \end{itemize}
 */

template<typename T, typename T_Storage>
gEigenData< std::complex<T>, std::complex<T> > eig(const gMatrix<std::complex<T>, T_Storage>& data);

/**
 * Return Eigenvectors and -values of a Matrix of {\bf complex<T>} into the
 * supplied argument matrices of {\bf complex<T> }.
 * @param data [N x N] Matrix of {\em complex<T>} containing data vectors.
 * @param evalues [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvectors as columns on exit. (Ignored on entry).
 * @param evectors [N x N] Matrix of {\em complex<T> } containing the
 *     eigenvalues on the diagonal on exit. (Ignored on entry).
 */
template<typename T, typename T_Storage1, typename T_Storage2, typename T_Storage3>
void eig(const gMatrix<std::complex<T>, T_Storage1>& data,
         gMatrix<std::complex<T>, T_Storage2>& evectors,
         gMatrix<std::complex<T>, T_Storage3>& evalues);
//@}

/**@name Singular Value Decomposition (SVD)
 */
//@{
/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf float } into the supplied argument matrices of
 * {\bf float }. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x N]   [N x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em float} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em float}.
 * @param S [M x N] Diagonal Matrix of {\em float}.
 * @param VT [N x N] Unitary Matrix of {\em float}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<float, T_Storage1>& A,
         gMatrix<float, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<float, T_Storage4>& VT);

/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf complex<float> } into the supplied argument matrices of
 * {\bf complex<float> } and {\bf float}. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x N]   [N x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em complex<float>} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em complex<float>}.
 * @param S [M x N] Diagonal Matrix of {\em float}.
 * @param VT [N x N] Unitary Matrix of {\em complex<float>}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<std::complex<float>, T_Storage1>& A,
         gMatrix<std::complex<float>, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<std::complex<float>, T_Storage4>& VT);

/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf double } into the supplied argument matrices of
 * {\bf double }. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x N]   [N x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em double} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em double}.
 * @param S [M x N] Diagonal Matrix of {\em double}.
 * @param VT [N x N] Unitary Matrix of {\em double}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<double, T_Storage1>& A,
         gMatrix<double, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<double, T_Storage4>& VT);

/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf complex<double> } into the supplied argument matrices of
 * {\bf complex<double> } and {\bf double}. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x N]   [N x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em complex<double>} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em complex<double>}.
 * @param S [M x N] Diagonal Matrix of {\em double}.
 * @param VT [N x N] Unitary Matrix of {\em complex<double>}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void svd(const gMatrix<std::complex<double>, T_Storage1>& A,
         gMatrix<std::complex<double>, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<std::complex<double>, T_Storage4>& VT);

//@}

/**@name Reduced Singular Value Decomposition (SVD)
 */
//@{
/**
 * Return the Reduced Singular Value Decomposition (SVD) of a matrix of
 * {\bf float } into the supplied argument matrices of
 * {\bf float }. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x M]   [M x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em float} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em float}.
 * @param S [M x M] Diagonal Matrix of {\em float}.
 * @param VT [M x N] Unitary Matrix of {\em float}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<float, T_Storage1>& A,
         gMatrix<float, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<float, T_Storage4>& VT);

/**
 * Return the Reduced Singular Value Decomposition (SVD) of a matrix of
 * {\bf complex<float> } into the supplied argument matrices of
 * {\bf complex<float> } and {\bf float}. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x M]   [M x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em complex<float>} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em complex<float>}.
 * @param S [M x M] Diagonal Matrix of {\em float}.
 * @param VT [M x N] Unitary Matrix of {\em complex<float>}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<std::complex<float>, T_Storage1>& A,
         gMatrix<std::complex<float>, T_Storage2>& U,
         gMatrix<float, T_Storage3>& S,
         gMatrix<std::complex<float>, T_Storage4>& VT);

/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf double } into the supplied argument matrices of
 * {\bf double }. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x M]   [M x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em double} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em double}.
 * @param S [M x M] Diagonal Matrix of {\em double}.
 * @param VT [M x N] Unitary Matrix of {\em double}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<double, T_Storage1>& A,
         gMatrix<double, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<double, T_Storage4>& VT);

/**
 * Return the Singular Value Decomposition (SVD) of a matrix of
 * {\bf complex<double> } into the supplied argument matrices of
 * {\bf complex<double> } and {\bf double}. \\
 * The input matrix is decomposed into three matrices as follows: \\
 * \begin{verbatim}
      A    =    U   *   S    *   VT

   [M x N]   [M x M] [M x M]   [M x N]

      U and VT are unitary and S is diagonal.
 \end{verbatim}
 * The matrix to be decomposed is passed as the first argument.
 * It is left untouched. The remaining arguments ( U, S and VT )
 * are resized to required size and the contents overwritten by the
 * function.
 * @param A [M x N] Matrix of {\em complex<double>} containing data vectors.
 * @param U [M x M] Unitary Matrix of {\em complex<double>}.
 * @param S [M x M] Diagonal Matrix of {\em double}.
 * @param VT [M x N] Unitary Matrix of {\em complex<double>}.
 */
template<typename T_Storage1, typename T_Storage2, typename T_Storage3, typename T_Storage4>
void r_svd(const gMatrix<std::complex<double>, T_Storage1>& A,
         gMatrix<std::complex<double>, T_Storage2>& U,
         gMatrix<double, T_Storage3>& S,
         gMatrix<std::complex<double>, T_Storage4>& VT);

//@}

GLINEAR_NAMESPACE_END

/*
 * Support functions.
 * All the following call the required Fortran function, setting up any support
 * variables and allocating workspace. This provides a clean interface to all the
 * Lapack routines without worrying about workspace allocation, passing
 * pointers or any other C-like nasties.
 */

// {s,c,d,z} potrf (cholesky factorization)
GLINEAR_API void o__potrf_dispatch(fortran::real* a,
                                   fortran::integer n,
                                   fortran::integer& info);

GLINEAR_API void o__potrf_dispatch(fortran::complex* a,
                                   fortran::integer n,
                                   fortran::integer& info);

GLINEAR_API void o__potrf_dispatch(fortran::doublereal* a,
                                   fortran::integer n,
                                   fortran::integer& info);

GLINEAR_API void o__potrf_dispatch(fortran::doublecomplex* a,
                                   fortran::integer n,
                                   fortran::integer& info);

GLINEAR_API void o__potrs_dispatch(fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::real* a,
                                   fortran::real* b,
                                   fortran::integer& info);

GLINEAR_API void o__potrs_dispatch(fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::complex* a,
                                   fortran::complex* b,
                                   fortran::integer& info);

GLINEAR_API void o__potrs_dispatch(fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::doublereal* a,
                                   fortran::doublereal* b,
                                   fortran::integer& info);

GLINEAR_API void o__potrs_dispatch(fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::doublecomplex* a,
                                   fortran::doublecomplex* b,
                                   fortran::integer& info);


// {s,c,d,z} getrf (LU factorization)
GLINEAR_API void o__getrf_dispatch(fortran::real* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getrf_dispatch(fortran::complex* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getrf_dispatch(fortran::doublereal* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getrf_dispatch(fortran::doublecomplex* a,
                                   fortran::integer m,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

// {s,c,d,z} laswp (Apply Permutation)
GLINEAR_API fortran::integer o__laswp_dispatch(fortran::real* a,
                                               fortran::integer m,
                                               fortran::integer n,
                                               fortran::integer* ipiv);

GLINEAR_API fortran::integer o__laswp_dispatch(fortran::complex* a,
                                               fortran::integer m,
                                               fortran::integer n,
                                               fortran::integer* ipiv);

GLINEAR_API fortran::integer o__laswp_dispatch(fortran::doublereal* a,
                                               fortran::integer m,
                                               fortran::integer n,
                                               fortran::integer* ipiv);

GLINEAR_API fortran::integer o__laswp_dispatch(fortran::doublecomplex* a,
                                               fortran::integer m,
                                               fortran::integer n,
                                               fortran::integer* ipiv);
// {s,c,d,z} getrs (LU Solve)
GLINEAR_API void o__getrs_dispatch(char trans,
                                   fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::real* a,
                                   fortran::integer* ipiv,
                                   fortran::real* b,
                                   fortran::integer& info);

GLINEAR_API void o__getrs_dispatch(char trans,
                                   fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::complex* a,
                                   fortran::integer* ipiv,
                                   fortran::complex* b,
                                   fortran::integer& info);

GLINEAR_API void o__getrs_dispatch(char trans,
                                   fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::doublereal* a,
                                   fortran::integer* ipiv,
                                   fortran::doublereal* b,
                                   fortran::integer& info);

GLINEAR_API void o__getrs_dispatch(char trans,
                                   fortran::integer n,
                                   fortran::integer nrhs,
                                   fortran::doublecomplex* a,
                                   fortran::integer* ipiv,
                                   fortran::doublecomplex* b,
                                   fortran::integer& info);


// {s,c,d,z} getri (matrix inversion)
GLINEAR_API void o__getri_dispatch(fortran::real* a,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getri_dispatch(fortran::complex* a,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getri_dispatch(fortran::doublereal* a,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

GLINEAR_API void o__getri_dispatch(fortran::doublecomplex* a,
                                   fortran::integer n,
                                   fortran::integer* ipiv,
                                   fortran::integer& info);

// {c,z} geevx (eigenanalysis expert driver)
GLINEAR_API fortran::integer o__geevx_dispatch(fortran::complex* a,
           fortran::integer n,
           fortran::complex* vr,
           fortran::complex* w);

GLINEAR_API fortran::integer o__geevx_dispatch(fortran::doublecomplex* a,
           fortran::integer n,
           fortran::doublecomplex* vr,
           fortran::doublecomplex* w);

// {s,c,d,z} gesvd (singular value decomposition)
GLINEAR_API fortran::integer o__gesvd_dispatch(fortran::real* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::real* u,
           fortran::real* vt);

GLINEAR_API fortran::integer o__gesvd_dispatch(fortran::complex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::complex* u,
           fortran::complex* vt);

GLINEAR_API fortran::integer o__gesvd_dispatch(fortran::doublereal* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublereal* u,
           fortran::doublereal* vt);

GLINEAR_API fortran::integer o__gesvd_dispatch(fortran::doublecomplex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublecomplex* u,
           fortran::doublecomplex* vt);

// reduced {s,c,d,z} gesvd (singular value decomposition)
GLINEAR_API fortran::integer o__reduced_gesvd_dispatch(fortran::real* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::real* u,
           fortran::real* vt);

GLINEAR_API fortran::integer o__reduced_gesvd_dispatch(fortran::complex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::real* s,
           fortran::complex* u,
           fortran::complex* vt);

GLINEAR_API fortran::integer o__reduced_gesvd_dispatch(fortran::doublereal* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublereal* u,
           fortran::doublereal* vt);

GLINEAR_API fortran::integer o__reduced_gesvd_dispatch(fortran::doublecomplex* a,
           fortran::integer m,
           fortran::integer n,
           fortran::doublereal* s,
           fortran::doublecomplex* u,
           fortran::doublecomplex* vt);

/*
  Support functions for pivot to permutation convertion
 */
void showipiv(const fortran::fortran_array<fortran::integer>& ipiv,
              int sz);

/**
 * transforms lapack ipiv to row permutation vector
 * both arrays MUST be presized
 */
void ipiv2perm(const fortran::fortran_array<fortran::integer>& ipiv,
               gLinear::gRowVector<int>& perm);

/**
 * transforms row permutation vector to lapack ipiv
 * both arrays must be presized
 */
void perm2ipiv(const gLinear::gRowVector<int>& perm,
               fortran::fortran_array<fortran::integer>& ipiv);


#include "gLapackWrap.tcc"

#endif // GLAPACK_WRAP_H
