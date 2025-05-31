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

// standard headers
#include <iostream>  // cout, endl
#include <cmath>  // ceil
#include <new>  // placement new


GLINEAR_NAMESPACE_BEGIN

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(void) :
  T_Base() {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(T_Index _rows,
               T_Index _cols) :
  T_Base(_rows, _cols) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(const T_Matrix& mat) :
  T_Base(mat) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(const gRowVector<T_Value, DENSE>& vec) :
  T_Base( 1, vec.size(), vec.getStor(), vec.getStartAddr() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(const gColVector<T_Value, DENSE>& vec) :
  T_Base( vec.size(), 1, vec.getStor(), vec.getStartAddr() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::gMatrixCore(T_ConstPointer addr,
                                                         T_Index rows, T_Index cols) :
  T_Base(addr, rows, cols, matrix_no_copy_tag()) {
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::configure(T_Index new_nRows,
                  T_Index new_nCols,
                  const T_StorageElem* new_stor,
                  T_ConstPointer new_startAddr) {
  T_Base::configure(new_nRows, new_nCols, new_stor, new_startAddr);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Pointer
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::mapAddr(T_Index _row,
                 T_Index _col) const {
  return const_cast<T_Pointer>( this->getStartAddr() + _row * this->cols() + _col );
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialInfo(int i) const {
  std::cout << "gMatrixCore< gMatrix<T_Value, ROW_DENSE> >" << std::endl;
  T_Base::specialInfo(i);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Transpose
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialTranspose(void) const {
  T_Transpose tmat = gLinearNoInit();
  tmat.configure( this->cols(), this->rows(), this->getStor(), this->getStartAddr() );
  return tmat;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Row
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialRow(T_Index _row) const {
  return T_Row( this->cols(), this->getStor(), mapAddr(_row, 0) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Col
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialCol(T_Index _col) const {
  return T_Col( this->rows(), this->getStor(), mapAddr(0, _col), this->cols() );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialRowBegin(void) const {
  return T_RowIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialRowEnd(void) const {
  return T_RowIterator( this->rows(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialColBegin(void) const {
  return T_ColIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialColEnd(void) const {
  return T_ColIterator( this->cols(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Reference
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialElem(T_Index _row,
               T_Index _col) const {
  return *mapAddr(_row, _col);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Rowsub
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialExtractRows(T_Index fromRow,
                T_Index toRow) const {
  T_Rowsub tmp = gLinearNoInit();
  tmp.configure( toRow - fromRow + 1, this->cols(), this->getStor(), mapAddr(fromRow, 0) );
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Colsub
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialExtractCols(T_Index fromCol,
                T_Index toCol) const {
  T_Colsub tmp = gLinearNoInit();
  tmp.configure( this->rows(), toCol - fromCol + 1, this->getStor(), mapAddr(0, fromCol), this->cols(), 1 );
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::T_Slice
gMatrixCore< gMatrix<TF_Value, ROW_DENSE> >::specialSlice(const gIndexRange& _row,
                const gIndexRange& _col) const {

  const T_Index new_nRows = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + _row.getEnd() - _row.getStart() ) /
                   static_cast<double>( _row.getStride() ) ) );
  const T_Index new_nCols = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + _col.getEnd() - _col.getStart() ) /
                   static_cast<double>( _col.getStride() ) ) );
  T_Slice tmp = gLinearNoInit();
  tmp.configure( new_nRows, new_nCols, this->getStor(),
     mapAddr( _row.getStart(), _col.getStart() ),
     this->cols() * _row.getStride(), _col.getStride() );
  return tmp;
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(void) :
  T_Base() {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(T_Index _rows,
               T_Index _cols) :
  T_Base(_rows, _cols) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(const T_Matrix& mat) :
  T_Base(mat) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(const gRowVector<T_Value, DENSE>& vec) :
  T_Base( 1, vec.size(), vec.getStor(), vec.getStartAddr() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(const gColVector<T_Value, DENSE>& vec) :
  T_Base( vec.size(), 1, vec.getStor(), vec.getStartAddr() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::gMatrixCore(T_ConstPointer addr,
                                                         T_Index rows, T_Index cols) :
  T_Base(addr, rows, cols, matrix_no_copy_tag()) {
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::configure(T_Index new_nRows,
                  T_Index new_nCols,
                  const T_StorageElem* new_stor,
                  T_ConstPointer new_startAddr) {
  T_Base::configure(new_nRows, new_nCols, new_stor, new_startAddr);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Pointer
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::mapAddr(T_Index _row,
                 T_Index _col) const {
  return const_cast<T_Pointer>(this->getStartAddr() + _col * this->rows() + _row);
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialInfo(int i) const {
  std::cout << "gMatrixCore< gMatrix<T_Value, COL_DENSE> >" << std::endl;
  T_Base::specialInfo(i);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Transpose
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialTranspose(void) const {
  T_Transpose tmat = gLinearNoInit();
  tmat.configure( this->cols(), this->rows(), this->getStor(), this->getStartAddr() );
  return tmat;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Row
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialRow(T_Index _row) const {
  return T_Row( this->cols(), this->getStor(), mapAddr(_row, 0), this->rows() );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Col
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialCol(T_Index _col) const {
  return T_Col( this->rows(), this->getStor(), mapAddr(0, _col) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialRowBegin(void) const {
  return T_RowIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialRowEnd(void) const {
  return T_RowIterator( this->rows(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialColBegin(void) const {
  return T_ColIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialColEnd(void) const {
  return T_ColIterator( this->cols(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Reference
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialElem(T_Index _row,
               T_Index _col) const {
  return *mapAddr(_row, _col);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Rowsub
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialExtractRows(T_Index fromRow,
                T_Index toRow) const {
  T_Rowsub tmp = gLinearNoInit();
  tmp.configure( toRow - fromRow + 1, this->cols(), this->getStor(), mapAddr(fromRow, 0), 1, this->rows() );
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Colsub
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialExtractCols(T_Index fromCol,
                T_Index toCol) const {
  T_Colsub tmp = gLinearNoInit();
  tmp.configure( this->rows(), toCol - fromCol + 1, this->getStor(), mapAddr(0, fromCol) );
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::T_Slice
gMatrixCore< gMatrix<TF_Value, COL_DENSE> >::specialSlice(const gIndexRange& _row,
                const gIndexRange& _col) const {
  const T_Index new_nRows = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + _row.getEnd() - _row.getStart() ) /
                   static_cast<double>( _row.getStride() ) ) );
  const T_Index new_nCols = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + _col.getEnd() - _col.getStart() ) /
                   static_cast<double>( _col.getStride() ) ) );
  T_Slice tmp = gLinearNoInit();
  tmp.configure( new_nRows, new_nCols, this->getStor(),
     mapAddr( _row.getStart(), _col.getStart() ),
     _row.getStride(), this->rows() * _col.getStride() );
  return tmp;
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(void) :
  T_Base() {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(T_Index rows,
                                                     T_Index cols) :
  T_Base(rows, cols),
  rowStride(cols),
  colStride(1) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(const T_Matrix& mat) :
  T_Base(mat),
  rowStride(mat.rowStride),
  colStride(mat.colStride) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(T_Index new_nRows,
                                                     T_Index new_nCols,
                                                     const T_StorageElem* new_stor,
                                                     T_ConstPointer new_startAddr,
                                                     T_Index new_rowStride,
                                                     T_Index new_colStride) :
  T_Base(new_nRows, new_nCols, new_stor, new_startAddr),
  rowStride(new_rowStride),
  colStride(new_colStride) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(const gRowVector<T_Value, DENSE>& vec) :
  T_Base(1, vec.size(), vec.getStor(), vec.getStartAddr(), vec.size(), 1) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(const gRowVector<T_Value, SLICE>& vec) :
  T_Base( 1, vec.size(), vec.getStor(), vec.getStartAddr(), vec.size(), vec.getStride * vec.size() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(const gColVector<T_Value, DENSE>& vec) :
  T_Base( vec.size(), 1, vec.getStor(), vec.getStartAddr(), 1, vec.size() ) {
}

template<typename TF_Value>
gMatrixCore< gMatrix<TF_Value, SLICE> >::gMatrixCore(const gColVector<T_Value, SLICE>& vec) :
  T_Base( vec.size(), 1, vec.getStor(), vec.getStartAddr(), vec.getStride() * vec.size(), vec.size() ) {
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, SLICE> >::configure(T_Index new_nRows,
                                                        T_Index new_nCols,
                                                        const T_StorageElem* new_stor,
                                                        T_ConstPointer new_startAddr,
                                                        T_Index new_rowStride,
                                                        T_Index new_colStride) {
  this->~gMatrixCore();
  new(this) gMatrixCore(new_nRows, new_nCols,
                        new_stor, new_startAddr,
                        new_rowStride, new_colStride);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Pointer
gMatrixCore< gMatrix<TF_Value, SLICE> >::mapAddr(T_Index _row,
             T_Index _col) const {
  return const_cast<T_Pointer>(this->getStartAddr() + _row * rowStride + _col * colStride);
}

template<typename TF_Value>
void gMatrixCore< gMatrix<TF_Value, SLICE> >::specialInfo(int i) const {
  std::cout << "gMatrixCore< gMatrix<T_Value, SLICE> >" << std::endl;
  T_Base::specialInfo(i);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Transpose
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialTranspose(void) const {
  T_Transpose tmp = gLinearNoInit();
  tmp.configure(this->cols(), this->rows(), this->getStor(), this->getStartAddr(), colStride, rowStride);
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Row
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialRow(T_Index _row) const {
  return T_Row(this->cols(), this->getStor(), mapAddr(_row, 0), colStride);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Col
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialCol(T_Index _col) const {
  return T_Col(this->rows(), this->getStor(), mapAddr(0, _col), rowStride);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialRowBegin(void) const {
  return T_RowIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_RowIterator
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialRowEnd(void) const {
  return T_RowIterator( this->rows(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialColBegin(void) const {
  return T_ColIterator( 0, reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_ColIterator
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialColEnd(void) const {
  return T_ColIterator( this->cols(), reinterpret_cast<const T_Matrix*>(this) );
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Reference
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialElem(T_Index _row,
                 T_Index _col) const {
  return *mapAddr(_row, _col);
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Rowsub
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialExtractRows(T_Index fromRow,
                  T_Index toRow) const {
  T_Rowsub tmp = gLinearNoInit();
  tmp.configure(toRow - fromRow + 1, this->cols(), this->getStor(), mapAddr(fromRow, 0), rowStride, colStride);
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Colsub
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialExtractCols(T_Index fromCol,
                  T_Index toCol) const {
  T_Colsub tmp = gLinearNoInit();
  tmp.configure(this->rows(), toCol - fromCol + 1, this->getStor(), mapAddr(0, fromCol), rowStride, colStride);
  return tmp;
}

template<typename TF_Value>
typename gMatrixCore< gMatrix<TF_Value, SLICE> >::T_Slice
gMatrixCore< gMatrix<TF_Value, SLICE> >::specialSlice(const gIndexRange& row,
                  const gIndexRange& col) const {
  const T_Index new_nRows = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + row.getEnd() - row.getStart() ) /
                   static_cast<double>( row.getStride() ) ) );
  const T_Index new_nCols = static_cast<T_Index>( std::ceil( static_cast<double>( 1 + col.getEnd() - col.getStart() ) /
                   static_cast<double>( col.getStride() ) ) );
  T_Slice tmp = gLinearNoInit();
  tmp.configure( new_nRows, new_nCols,
                 this->getStor(), mapAddr( row.getStart(), col.getStart() ),
                 rowStride * row.getStride(), colStride * col.getStride() );
  return tmp;
}

GLINEAR_NAMESPACE_END
