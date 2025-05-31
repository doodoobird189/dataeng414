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

#ifndef GMATRIX_BASE_H
#define GMATRIX_BASE_H

// glinear headers
#include "gIndexRange.h"
#include "gMatrixTraits.h"
#include "gLinearMacros.h"

// standard headers
#include <iosfwd>  // ostream, istream
#include <string>  // string

GLINEAR_NAMESPACE_BEGIN

/**
 * Generic Matrix base class.
 */
template<typename TF_Matrix>
class gMatrixBase {
public:
  /**
   * Type of this matrix.
   */
  typedef TF_Matrix T_Matrix;
  /**
   * Type of matrix elements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Value T_Value;
  /**
   * Type of storage the matrix implements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Storage T_Storage;
  /**
   * Type of const pointer to matrix elements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstPointer T_ConstPointer;
  /**
   * Type of pointer to matrix elements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Pointer T_Pointer;
  /**
   * Type of const reference to matrix elements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstReference T_ConstReference;
  /**
   * Type of reference to matrix elements.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Reference T_Reference;
  /**
   * Type of const matrix slice.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstSlice T_ConstSlice;
  /**
   * Type of matrix slice.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Slice T_Slice;
  /**
   * Type of const row.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstRow T_ConstRow;
  /**
   * Type of row.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Row T_Row;
  /**
   * Type of const column.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstCol T_ConstCol;
  /**
   * Type of column.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Col T_Col;
  /**
   * Const row iterator.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstRowIterator T_ConstRowIterator;
  /**
   * Row iterator.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_RowIterator T_RowIterator;
  /**
   * Const column iterator.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstColIterator T_ConstColIterator;
  /**
   * Column iterator.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ColIterator T_ColIterator;
  /**
   * Type of contiguous const row subrange.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstRowsub T_ConstRowsub;
  /**
   * Type of contiguous row subrange.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Rowsub T_Rowsub;
  /**
   * Type of contiguous const column subrange.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstColsub T_ConstColsub;
  /**
   * Type of contiguous column subrange.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Colsub T_Colsub;
  /**
   * Index type
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Index T_Index;
  /**
   * Type of const transpose matrix.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_ConstTranspose T_ConstTranspose;
  /**
   * Type of transpose matrix.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Transpose T_Transpose;
  /**
   * Type of dense matrix.
   */
  typedef typename gMatrixTraits<T_Matrix>::T_Dense T_Dense;

  /**
   * ID String
   */
  static const char className[];
  /**
   * Info dump starting at most derived
   */
  void info(int i = 0) const;
  /**
   * Return number of rows in matrix.
   * @return Number of rows.
   */
  T_Index rows(void) const;
  /**
   * Return number of columns in matrix.
   * @return Number of columns.
   */
  T_Index cols(void) const;
  /**
   * Resize the matrix. \\
   * Resize the matrix to [{\em rows} x {\em cols}]. \\
   * If CHK_IDX is defined, the resize is only performed if this is
   * the only object using the storage class (i.e. reference count
   * equal to 1). \\
   * This call returns without doing anything if the
   * requested size equals the current size. \\
   * {\bf Note:} All the data is lost as a new RefCntStorage object
   * is allocated. This means that this matrix is first
   * {\em depeened}, and then resized.
   * This leaves all previous aliases untouched.
   * This is in effect the same as calling the
   * Matrix(int, int) constructor.
   * @param rows New number of rows.
   * @param cols New number of columns.
   * @return
   * \begin{itemize}
   * \item {\bf true - } Resize successfully executed
   * \item {\bf false - } Resize returned without modifying this object.
   * \end{itemize}
   */
  bool resize(T_Index new_rows,
              T_Index new_cols);
protected:
  /**
   * General Matrix - Matrix assigment operator. \\
   * Resize this matrix and initialize with the values of the argument
   * matrix {\em mat}.
   * A {\bf deep} copy is performed from {\em mat} to this Matrix.
   * @param mat Source Matrix.
   * @return *this.
   */
  template<typename T_SrcMatrix>
  T_Matrix& assignTo(const gMatrixBase<T_SrcMatrix>& mat);

public:

  /**
   * Scalar assignment. \\
   * Assign the argument value {\em val} to all of the elements in
   * the Matrix.
   * @param scalar Value to assign to all elements.
   * @return *this.
   */
  T_Matrix& assignToAll(T_ConstReference scalar);

  /**
   * Return const transpose of this matrix. \\
   * @return Const transpose of this matrix.
   */
  T_ConstTranspose transpose(void) const;
  /**
   * Return transpose of this Matrix. \\
   * @return Transpose of this Matrix.
   */
  T_Transpose transpose(void);
  /**
   * Return const gRowVector representing a single row of the Matrix.
   * Optional (CHK_IDX) index bound checking is performed.
   * @param row Row index (0 to nRows-1).
   * @return gRowVector representing the {\em row}'th row.
   */
  T_ConstRow row(T_Index row) const;
  /**
   * Return gRowVector representing a single row of the Matrix. \\
   * Optional (CHK_IDX) index bound checking is performed.
   * @param row Row index (0 to nRows-1).
   * @return gRowVector<T, DENSE> representing the {\em row}'th row.
   */
  T_Row row(T_Index row);
  /**
   * Return const gColVector representing a single column of the Matrix.
   * Optional (CHK_IDX) index bound checking is performed.
   * @param col Column index (0 to nCols-1).
   * @return gColVector representing the {\em col}'th column.
   */
  T_ConstCol col(T_Index col) const;
  /**
   * Return gColVector representing a single column of the Matrix. \\
   * Optional (CHK_IDX) index bound checking is performed.
   * @param col Column index (0 to nCols-1).
   * @return ColVector<T, SLICE> representing the {\em col}'th column..
   */
  T_Col col(T_Index col);
  /**
   * Return const row begin iterator. \\
   * Dereferncing this iterator returns the first row.
   * The iterator may be moved to point to other rows through the
   * increment and decrement operators.
   */
  T_ConstRowIterator rowBegin(void) const;
  /**
   * Return row begin iterator. \\
   * Dereferncing this iterator returns the first row.
   * The iterator may be moved to point to other rows through the
   * increment and decrement operators.
   */
  T_RowIterator rowBegin(void);
  /**
   * Return const row past-the-end iterator. \\
   * Dereferncing this iterator is illegal.
   * The iterator may be moved to point to other rows through the
   * increment and decrement operators.
   */
  T_ConstRowIterator rowEnd(void) const;
  /**
   * Return row past-the-end iterator. \\
   * Dereferncing this iterator is illegal.
   * The iterator may be moved to point to other rows through the
   * increment and decrement operators.
   */
  T_RowIterator rowEnd(void);
  /**
   * Return const column begin iterator. \\
   * Dereferncing this iterator returns the first column.
   * The iterator may be moved to point to other columns through the
   * increment and decrement operators.
   */
  T_ConstColIterator colBegin(void) const;
  /**
   * Return column begin iterator. \\
   * Dereferncing this iterator returns the first column.
   * The iterator may be moved to point to other columns through the
   * increment and decrement operators.
   */
  T_ColIterator colBegin(void);
  /**
   * Return const column past-the-end iterator. \\
   * Dereferncing this iterator is illegal.
   * The iterator may be moved to point to other columns through the
   * increment and decrement operators.
   */
  T_ConstColIterator colEnd(void) const;
  /**
   * Return column past-the-end iterator. \\
   * Dereferncing this iterator is illegal.
   * The iterator may be moved to point to other columns through the
   * increment and decrement operators.
   */
  T_ColIterator colEnd(void);
  /**
   * Alias of const T_Row row(T_Index row) const.
   */
  T_ConstRow operator[](T_Index row) const;
  /**
   * Alias of T_Row row(T_Index row).
   */
  T_Row operator[](T_Index row);
  /**
   * Constant element access operator. \\
   * Return const reference to the matrix element in row {\em row} and
   * column {\em col}. \\
   * @param Element row.
   * @param Element column.
   * @return Const reference to element ({\em row}, {\em col})
   * as non-modifyable lvalue.
   */
  T_ConstReference operator()(T_Index row, T_Index col) const;
  /**
   * Element access operator. \\
   * Return reference to the matrix element in row {\em row} and column
   * {\em col}. \\
   * @param row Element row.
   * @param col Element column.
   * @return Reference to element ({\em row}, {\em col}) as modifyable lvalue.
   */
  T_Reference operator()(T_Index row, T_Index col);

  /**
   * Extract a contiguous const subrange of rows as a
   * const T_Rowsub matrix. \\
   * {\bf Note: } The extraction returns a {\bf shallow} copy of the
   * selected elements. This means that any changes made to the
   * extracted matrix would reflect in the original matrix.
   * @param fromRow Start row index.
   * @param toRow End row index.
   * @return New matrix containing a subrange of rows \\
   * from this matrix.
   */
  T_ConstRowsub extractRows(T_Index fromRow,
                            T_Index toRow) const;
  /**
   * Extract a contiguous subrange of rows as a T_Rowsub matrix. \\
   * {\bf Note: } The extraction returns a {\bf shallow} copy of the
   * selected elements. This means that any changes made to the
   * extracted matrix would reflect in the original matrix.
   * @param fromRow Start row index.
   * @param toRow End row index.
   * @return New matrix containing a subrange of rows \\
   * from this matrix.
   */
  T_Rowsub extractRows(T_Index fromRow,
                       T_Index toRow);
  /**
   * Extract a contiguous const subrange of columns as a
   * const T_Colsub Matrix. \\
   * {\bf Note: } The extraction returns a {\bf shallow} copy of the
   * selected elements. This means that any changes made to the
   * extracted matrix would reflect in the original matrix.
   * @param fromCol Start row index.
   * @param toCol End row index.
   * @return New slice Matrix containing a subrange of columns \\
   * from this Matrix.
   *
   */
  T_ConstColsub extractCols(T_Index fromCol,
                            T_Index toCol) const;
  /**
   * Extract a contiguous subrange of columns as a T_Colsub Matrix. \\
   * {\bf Note: } The extraction returns a {\bf shallow} copy of the
   * selected elements. This means that any changes made to the
   * extracted matrix would reflect in the original matrix.
   * @param fromCol Start row index.
   * @param toCol End row index.
   * @return New slice Matrix containing a subrange of columns \\
   * from this Matrix.
   *
   */
  T_Colsub extractCols(T_Index fromCol,
                       T_Index toCol);
  /**
   * Return complete matrix as constant slice. \\
   */
  T_ConstSlice slice(void) const;
  /**
   * Return complete matrix as slice. \\
   */
  T_Slice slice(void);
  /**
   * Return subrange slice of Matrix. \\
   * The elements specified by the {\em rows} and {\em cols} arguments
   * determine the slice of the matrix that is returned. \\
   * {\bf Note: } The returned slice is a {\bf shallow} copy of this
   * matrix.
   * This implies that any changes made to the elements of this slice
   * affects this original matrix.
   * @param rows Range (start, end, stride) of rows to return.
   * @param cols Range (start, end, stride) of columns to return.
   * @return Slice representing the selected elements.
   */
  T_ConstSlice slice(const gIndexRange& rows,
                     const gIndexRange& cols) const;
  /**
   * Return subrange slice of Matrix. \\
   * The elements specified by the {\em rows} and {\em cols} arguments
   * determine the slice of the matrix that is returned. \\
   * {\bf Note: } The returned slice is a {\bf shallow} copy of this
   * matrix.
   * This implies that any changes made to the elements of this slice
   * affects this original matrix.
   * @param rows Range (start, end, stride) of rows to return.
   * @param cols Range (start, end, stride) of columns to return.
   * @return Slice representing the selected elements.
   */
  T_Slice slice(const gIndexRange& rows,
                const gIndexRange& cols);

  /**
   * Return the number of references (aliases) existing on this object.
   */
  T_Index getRefs(void) const;
  /**
   * Deepen.
   */
  void deepen(void);
  /**
   * Const alias.
   */
  void alias(const T_Matrix& mat) const;
  /**
   * Alias.
   */
  void alias(T_Matrix& mat);
  /**
   * Print the contents to stdout.
   */
  void print(int depth = 0) const;
  /**
   * Binary Output
   * Write the length of the vector, followed by the contets
   * to a stream (binary output).
   */
  std::ostream& write(std::ostream& os) const;
  /**
   * Binary Input
   * Read the length of the vector, followed by the contets
   * from a stream (binary input).
   */
  std::istream& read(std::istream& is);
private: // support
  void chkidxIndex(std::string funcName, std::string infostr, T_Index ubound, T_Index value) const;
  void chkidxIndex(std::string funcName, T_Index urbound, T_Index ucbound, T_Index row, T_Index col) const;
  void chkidxSubrange(std::string funcName, std::string infostr, T_Index from, T_Index to, T_Index ubound) const;
  void chkidxSlice(std::string funcName, const gIndexRange& rows, const gIndexRange& cols) const;
  void chkidxAlias(std::string funcName, std::string infostr) const;
};

GLINEAR_NAMESPACE_END

#include "gMatrixBase.tcc"

#endif // GMATRIX_BASE_H
