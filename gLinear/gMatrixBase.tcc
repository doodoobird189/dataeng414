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

// glinear headers
#include "gLinearMacros.h"

// standard headers
#include <iostream>  // cout, ostream, istream
#include <string>  // string


GLINEAR_NAMESPACE_BEGIN

template<typename T>
const char gMatrixBase<T>::className[] = "gMatrixBase<T>";


// dump info
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::info(int i) const {
  GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialInfo(i);
}

// return number of rows.
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Index
gMatrixBase<T_Matrix>::rows(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialRows();
}

// return number of columns.
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Index
gMatrixBase<T_Matrix>::cols(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialCols();
}

// resize matrix. decr ref cnt and delete if ref cnt zero.
// set nRows. set nCols. allocate new storage (rows*cols). inc ref cnt.
template<typename T_Matrix>
inline bool
gMatrixBase<T_Matrix>::resize(T_Index new_rows,
            T_Index new_cols) {
  if ( (rows() == new_rows) && (cols() == new_cols) ) {
    return false;
  }
#ifdef CHK_IDX
  // check for aliases
  chkidxAlias("resize(T_Index new_rows, T_Index new_cols)", "resize");
#endif // CHK_IDX
  GLINEAR_STATIC_VIRTUAL(T_Matrix) specialResize(new_rows, new_cols);
  return true;
}

// general matrix assignment.
// resize and deep copy from arg. vector to this.
template<typename T_Matrix> template<typename T_SrcMatrix>
inline typename gMatrixBase<T_Matrix>::T_Matrix&
gMatrixBase<T_Matrix>::assignTo(const gMatrixBase<T_SrcMatrix>& mat) {
  const T_Index nr = mat.rows();
  const T_Index nc = mat.cols();
  // attempt resize
  resize(nr, nc);
  // copy contents
  for (T_Index rindx = 0; rindx < nr; rindx++) {
    for (T_Index cindx = 0; cindx < nc; cindx++) {
      // most general and saves construction.
      operator()(rindx, cindx) = static_cast<T_Value>(mat(rindx, cindx));
    }
  }
  return reinterpret_cast<T_Matrix&>(*this); // downcast to md matrix class
}

// scalar assignment
template<typename T_Matrix>
typename gMatrixBase<T_Matrix>::T_Matrix&
gMatrixBase<T_Matrix>::assignToAll(T_ConstReference scalar) {
  const T_Index nr = rows();
  const T_Index nc = cols();
  // copy contents
  for (T_Index rindx = 0; rindx < nr; rindx++) {
    for (T_Index cindx = 0; cindx < nc; cindx++) {
      // most general and saves construction.
      operator()(rindx, cindx) = scalar;
    }
  }
  return reinterpret_cast<T_Matrix&>(*this); // downcast to md matrix class
}

// const transpose
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstTranspose
gMatrixBase<T_Matrix>::transpose(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialTranspose();
}

// transpose
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Transpose
gMatrixBase<T_Matrix>::transpose(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialTranspose();
}

// return const row
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstRow
gMatrixBase<T_Matrix>::row(T_Index _row) const {
#ifdef CHK_IDX
  chkidxIndex("row(T_Index row) const", "Invalid row index", rows() - 1, _row);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialRow(_row);
}

// return row
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Row
gMatrixBase<T_Matrix>::row(T_Index _row) {
#ifdef CHK_IDX
  chkidxIndex("row(T_Index row)", "Invalid row index", rows() - 1, _row);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialRow(_row);
}

// return const col
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstCol
gMatrixBase<T_Matrix>::col(T_Index _col) const {
#ifdef CHK_IDX
  chkidxIndex("col(T_Index col) const", "Invalid column index", cols() - 1, _col);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialCol(_col);
}

// return col
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Col
gMatrixBase<T_Matrix>::col(T_Index _col) {
#ifdef CHK_IDX
  chkidxIndex("col(T_Index col)", "Invalid column index", cols() - 1, _col);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialCol(_col);
}

// return const row begin iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstRowIterator
gMatrixBase<T_Matrix>::rowBegin(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialRowBegin();
}

// return row begin iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_RowIterator
gMatrixBase<T_Matrix>::rowBegin(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialRowBegin();
}

// return const row end iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstRowIterator
gMatrixBase<T_Matrix>::rowEnd(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialRowEnd();
}

// return row end iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_RowIterator
gMatrixBase<T_Matrix>::rowEnd(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialRowEnd();
}

// return const col begin iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstColIterator
gMatrixBase<T_Matrix>::colBegin(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialColBegin();
}

// return col begin iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ColIterator
gMatrixBase<T_Matrix>::colBegin(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialColBegin();
}

// return const col end iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstColIterator
gMatrixBase<T_Matrix>::colEnd(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialColEnd();
}

// return col end iterator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ColIterator
gMatrixBase<T_Matrix>::colEnd(void) {
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialColEnd();
}

// const row operator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstRow
gMatrixBase<T_Matrix>::operator[](T_Index r) const {
  return row(r);
}

// row operator
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Row
gMatrixBase<T_Matrix>::operator[](T_Index r) {
  return row(r);
}

// const element access. return const ref to element at pos (row, col).
// computed index.
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstReference
gMatrixBase<T_Matrix>::operator()(T_Index _row,
                                  T_Index _col) const {
#ifdef CHK_IDX
  chkidxIndex("operator()(T_Index row, T_Index col) const",
              rows() - 1, cols() - 1, _row, _col);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialElem(_row, _col);
}

// element access. return ref to element at pos (row, col). computed index.
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Reference
gMatrixBase<T_Matrix>::operator()(T_Index _row,
                                  T_Index _col) {
#ifdef CHK_IDX
  chkidxIndex("operator()(T_Index row, T_Index col)",
              rows() - 1, cols() - 1, _row, _col);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialElem(_row, _col);
}

// extract const rows into const T_Rowsub
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstRowsub
gMatrixBase<T_Matrix>::extractRows(T_Index fromRow,
                                   T_Index toRow) const {
#ifdef CHK_IDX
  chkidxSubrange("extractRows(T_Index fromRow, T_Index toRow) const", "",
                 fromRow, toRow, rows() - 1);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialExtractRows(fromRow, toRow);
}

// extract rows into T_Rowsub
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Rowsub
gMatrixBase<T_Matrix>::extractRows(T_Index fromRow,
                                   T_Index toRow) {
#ifdef CHK_IDX
  chkidxSubrange("extractRows(T_Index fromRow, T_Index toRow)", "",
                 fromRow, toRow, rows() - 1);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialExtractRows(fromRow, toRow);
}

// extract const cols into const T_Colsub
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstColsub
gMatrixBase<T_Matrix>::extractCols(T_Index fromCol,
                                   T_Index toCol) const {
#ifdef CHK_IDX
  chkidxSubrange("extractCols(T_Index fromCol, T_Index toCol) const", "",
                 fromCol, toCol, cols() - 1);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialExtractCols(fromCol, toCol);
}

// extract cols into T_Colsub
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Colsub
gMatrixBase<T_Matrix>::extractCols(T_Index fromCol,
                                   T_Index toCol) {
#ifdef CHK_IDX
  chkidxSubrange("extractCols(T_Index fromCol, T_Index toCol)", "",
                 fromCol, toCol, cols() - 1);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialExtractCols(fromCol, toCol);
}

// const full slice
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstSlice
gMatrixBase<T_Matrix>::slice(void) const {
  return slice( gIndexRange(0, rows() - 1), gIndexRange(0, cols() - 1) );
}

// full slice
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Slice
gMatrixBase<T_Matrix>::slice(void) {
  return slice( gIndexRange(0, rows() - 1), gIndexRange(0, cols() - 1) );
}

// const slice
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_ConstSlice
gMatrixBase<T_Matrix>::slice(const gIndexRange& _rows,
                             const gIndexRange& _cols) const {
#ifdef CHK_IDX
  chkidxSlice("slice(const gIndexRange& rows, const gIndexRange& cols) const",
              _rows, _cols);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialSlice(_rows, _cols);
}

// slice
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Slice
gMatrixBase<T_Matrix>::slice(const gIndexRange& _rows,
                             const gIndexRange& _cols) {
#ifdef CHK_IDX
  chkidxSlice("slice(const gIndexRange& rows, const gIndexRange& cols)",
              _rows, _cols);
#endif // CHK_IDX
  return GLINEAR_STATIC_VIRTUAL(T_Matrix) specialSlice(_rows, _cols);
}

// return number of refs
template<typename T_Matrix>
inline typename gMatrixBase<T_Matrix>::T_Index
gMatrixBase<T_Matrix>::getRefs(void) const {
  return GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialGetRefs();
}

// deepen
template<typename T_Matrix>
inline void gMatrixBase<T_Matrix>::deepen(void) {
  // return if this is the only view
  if ( getRefs() == 1 ) {
    return;
  }
  GLINEAR_STATIC_VIRTUAL(T_Matrix) specialDeepen();
}

// alias from const
template<typename T_Matrix>
inline void gMatrixBase<T_Matrix>::alias(const T_Matrix& mat) const {
  GLINEAR_STATIC_VIRTUAL(const T_Matrix) specialAlias(mat);
}

// alias
template<typename T_Matrix>
inline void gMatrixBase<T_Matrix>::alias(T_Matrix& mat) {
  GLINEAR_STATIC_VIRTUAL(T_Matrix) specialAlias(mat);
}

// ascii output
template<typename T_Matrix>
inline void
gMatrixBase<T_Matrix>::print(int depth) const {
  if (depth == 0) {
    return;
  } // if
  std::cout << reinterpret_cast<const T_Matrix&>(*this);
}

// binary output
template<typename T_Matrix>
inline std::ostream&
gMatrixBase<T_Matrix>::write(std::ostream& os) const {
  return global_write( os, reinterpret_cast<const T_Matrix&>(*this) );
}

// binary input
template<typename T_Matrix>
inline std::istream&
gMatrixBase<T_Matrix>::read(std::istream& is) {
  return global_read( is, reinterpret_cast<T_Matrix&>(*this) );
}

// 1D index bound checker
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::chkidxIndex(std::string funcName,
          std::string infostr,
          T_Index ubound,
          T_Index value) const {
  if (value > ubound) {
    gErrorIndex(className, funcName, infostr, ubound, value);
  }
}

// 2D index bound checker
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::chkidxIndex(std::string funcName,
                                        T_Index urbound,
          T_Index ucbound,
                                        T_Index _row,
          T_Index _col) const {
  if ( (_row > urbound) || (_col > ucbound) ) {
    gErrorIndex(className, funcName, urbound, ucbound, _row, _col);
  }
}

// subrange index checker
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::chkidxSubrange(std::string funcName,
             std::string infostr,
             T_Index from,
             T_Index to,
                                           T_Index ubound) const {
  if ( (from > to) || (to > ubound) ) {
    gErrorSubrange(className, funcName, infostr,
                   from, to, ubound);
  }
}

// slice index checker
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::chkidxSlice(std::string funcName,
                                        const gIndexRange& _row,
          const gIndexRange& _col) const {
  if ( ( _row.getStart() > _row.getEnd() ) || (_row.getEnd() > rows() - 1) ) {
    gErrorSlice( className, funcName, "Invalid row slice",
                 rows() - 1, _row.getStart(), _row.getEnd(), _row.getStride() );
  }
  if ( ( _col.getStart() > _col.getEnd() ) || (_col.getEnd() > cols() - 1) ) {
    gErrorSlice( className, funcName, "Invalid column slice",
                 cols() - 1, _col.getStart(), _col.getEnd(), _col.getStride() );
  }
}

// alias checker
template<typename T_Matrix>
void gMatrixBase<T_Matrix>::chkidxAlias(std::string funcName,
          std::string infostr) const {
  if (getRefs() > 1) {
    gErrorAliases(className, funcName, infostr);
  }
}


#if 0
template<typename T_Matrix>
RowVector<T, DENSE>*
gMatrixBase<T_Matrix>::arrayAddress(void) {
  // alloc storage for each RowVector
  // (each representing a row of the matrix)
  arrayAddress_storage = new RowVector<T, DENSE>[nRows];
  for (T_Index indx = 0; indx < nRows; indx++) {
    arrayAddress_storage[indx].alias( row(indx) );
  }
  return (RowVector<T, DENSE>*) arrayAddress_storage;
}

template<typename T_Matrix
const RowVector<T, DENSE>*
gMatrixBase<T_Matrix>::arrayAddress(void) const {
  // alloc storage for each RowVector
  // (each representing a row of the matrix)
  const RowVector<T, DENSE>* tmp =
    (const RowVector<T, DENSE>*) new RowVector<T, DENSE>[nRows];
  arrayAddress_storage = tmp;

  for (T_Index indx = 0; indx < nRows; indx++) {
    arrayAddress_storage[indx].alias( row(indx) );
  }
  return (const RowVector<T, DENSE>*) arrayAddress_storage;
}
#endif // 0

GLINEAR_NAMESPACE_END
