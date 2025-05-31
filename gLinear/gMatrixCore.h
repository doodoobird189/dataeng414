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

#ifndef GMATRIX_CORE_H
#define GMATRIX_CORE_H

// glinear headers
#include "gLinearMacros.h"
#include "gMatrixTraits.h"
#include "gMatrixMemBase.h"
#include "gMatrixIterator.h"
#include "gRowVector.h"
#include "gColVector.h"


GLINEAR_NAMESPACE_BEGIN

// forward decl. of matrix class
template<typename T_Value, typename T_Storage>
class gMatrix;

// primary template decl. of matrix core class
template<typename T_Matrix>
class gMatrixCore;

/**
 * Traits class for dense row-major matrix (ROW_DENSE)
 */
template<typename TF_Value>
class gMatrixTraits< gMatrix<TF_Value, ROW_DENSE> > {
public:
  /**
   * Type of matrix
   */
  typedef gMatrix<TF_Value, ROW_DENSE> T_Matrix;
  /**
   * Type of matrix elements
   */
  typedef TF_Value T_Value;
  /**
   * Type of const (immuting) pointer to elements
   */
  typedef const T_Value* T_ConstPointer;
  /**
   * Type of pointer to elements
   */
  typedef T_Value* T_Pointer;
  /**
   * Type of const (immuting) reference to elements
   */
  typedef const T_Value& T_ConstReference;
  /**
   * Type of reference to elements
   */
  typedef T_Value& T_Reference;
  /**
   * Index type
   */
  typedef typename gLinearTypes::T_Index T_Index;
  /**
   * Type of storage the implements
   */
  typedef ROW_DENSE T_Storage;
  /**
   * Const row
   */
  typedef const gRowVector<T_Value, DENSE> T_ConstRow;
  /**
   * Row
   */
  typedef gRowVector<T_Value, DENSE> T_Row;
  /**
   * Const column
   */
  typedef const gColVector<T_Value, SLICE> T_ConstCol;
  /**
   * Column
   */
  typedef gColVector<T_Value, SLICE> T_Col;
  /**
   * Const row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_ConstRowIterator;
  /**
   * Row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_RowIterator;
  /**
   * Const column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ConstColIterator;
  /**
   * Column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ColIterator;
  /**
   * Const contiguous row subrange
   */
  typedef const gMatrix<T_Value, ROW_DENSE> T_ConstRowsub;
  /**
   * Contiguous row subrange
   */
  typedef gMatrix<T_Value, ROW_DENSE> T_Rowsub;
  /**
   * Const contiguous column subrange
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstColsub;
  /**
   * Contiguous column subrange
   */
  typedef gMatrix<T_Value, SLICE> T_Colsub;
  /**
   * Const slice
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstSlice;
  /**
   * Slice
   */
  typedef gMatrix<T_Value, SLICE> T_Slice;
  /**
   * Dense type
   */
  typedef gMatrix<T_Value, ROW_DENSE> T_Dense;
  /**
   * Const transpose
   */
  typedef const gMatrix<T_Value, COL_DENSE> T_ConstTranspose;
  /**
   * Transpose
   */
  typedef gMatrix<T_Value, COL_DENSE> T_Transpose;
}; // gMatrixTraits


/**
 * Matrix core class spcialization for dense row-major matrices (ROW_DENSE)
 */
template<typename TF_Value>
class gMatrixCore< gMatrix<TF_Value, ROW_DENSE> > : public gMatrixMemBase< gMatrix<TF_Value, ROW_DENSE> > {
public:
  typedef gMatrix<TF_Value, ROW_DENSE> T_Matrix;
  typedef gMatrixMemBase<T_Matrix> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_Storage T_Storage;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_ConstReference T_ConstReference;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_ConstSlice T_ConstSlice;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_ConstRow T_ConstRow;
  typedef typename T_Base::T_Row T_Row;
  typedef typename T_Base::T_ConstCol T_ConstCol;
  typedef typename T_Base::T_Col T_Col;
  typedef typename T_Base::T_ConstRowIterator T_ConstRowIterator;
  typedef typename T_Base::T_RowIterator T_RowIterator;
  typedef typename T_Base::T_ConstColIterator T_ConstColIterator;
  typedef typename T_Base::T_ColIterator T_ColIterator;
  typedef typename T_Base::T_ConstRowsub T_ConstRowsub;
  typedef typename T_Base::T_Rowsub T_Rowsub;
  typedef typename T_Base::T_ConstColsub T_ConstColsub;
  typedef typename T_Base::T_Colsub T_Colsub;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;
  typedef typename T_Base::T_StorageElem T_StorageElem;
  typedef typename T_Base::T_Dense T_Dense;

public:
  // no init ctor
  inline gMatrixCore(gLinearNoInit);

  // default ctor
  inline gMatrixCore(void);

  // ctor specifying size
  inline gMatrixCore(T_Index rows,
                     T_Index cols);

   // copy ctor
  inline gMatrixCore(const T_Matrix& mat);

  // ctor from dense gRowVector
  inline gMatrixCore(const gRowVector<T_Value, DENSE>& vec);

  // ctor from dense gColVector
  inline gMatrixCore(const gColVector<T_Value, DENSE>& vec);

  /**
   * Initializing ctor 2: use existing space read-only.
   */
  gMatrixCore(T_ConstPointer new_startAddr, T_Index rows, T_Index cols);

  // configure with supplied parameters
  inline void configure(T_Index new_nRows,
                        T_Index new_nCols,
                        const T_StorageElem* new_stor,
                        T_ConstPointer new_startAddr);

  // map const address from matrix coords
  inline T_Pointer mapAddr(T_Index row,
                           T_Index col) const;

  // dump info
  inline void specialInfo(int i) const;

  // transpose
  inline T_Transpose specialTranspose(void) const;

  // return row
  inline T_Row specialRow(T_Index row) const;

  // return col
  inline T_Col specialCol(T_Index col) const;

  // return row begin iterator
  inline T_RowIterator specialRowBegin(void) const;

  // return row end iterator
  inline T_RowIterator specialRowEnd(void) const;

  // return col begin iterator
  inline T_ColIterator specialColBegin(void) const;

  // return col end iterator
  inline T_ColIterator specialColEnd(void) const;

  // return reference to element
  inline T_Reference specialElem(T_Index row,
                                 T_Index col) const;

  // return dense row subrange
  inline T_Rowsub specialExtractRows(T_Index fromRow,
                                     T_Index toRow) const;

  // return dense col subrange
  inline T_Colsub specialExtractCols(T_Index fromCol,
                                     T_Index toCol) const;

  // return slice
  inline T_Slice specialSlice(const gIndexRange& row,
                              const gIndexRange& col) const;

}; // gMatrixCore


/**
 * Traits class for dense column-major matrix (COL_DENSE)
 */
template<typename TF_Value>
class gMatrixTraits< gMatrix<TF_Value, COL_DENSE> > {
public:
  /**
   * Type of matrix.
   */
  typedef gMatrix<TF_Value, COL_DENSE> T_Matrix;
  /**
   * Type of matrix elements.
   */
  typedef TF_Value T_Value;
  /**
   * Type of const (immuting) pointer to elements.
   */
  typedef const T_Value* T_ConstPointer;
  /**
   * Type of pointer to elements.
   */
  typedef T_Value* T_Pointer;
  /**
   * Type of const (immuting) reference to elements.
   */
  typedef const T_Value& T_ConstReference;
  /**
   * Type of reference to elements.
   */
  typedef T_Value& T_Reference;
  /**
   * Index type.
   */
  typedef typename gLinearTypes::T_Index T_Index;
  /**
   * Type of storage the implements.
   */
  typedef COL_DENSE T_Storage;
  /**
   * Const row
   */
  typedef const gRowVector<T_Value, SLICE> T_ConstRow;
  /**
   * Row
   */
  typedef gRowVector<T_Value, SLICE> T_Row;
  /**
   * Const column
   */
  typedef const gColVector<T_Value, DENSE> T_ConstCol;
  /**
   * Column
   */
  typedef gColVector<T_Value, DENSE> T_Col;
  /**
   * Const row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_ConstRowIterator;
  /**
   * Row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_RowIterator;
  /**
   * Const column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ConstColIterator;
  /**
   * Column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ColIterator;
  /**
   * Const contiguous row subrange.
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstRowsub;
  /**
   * Contiguous row subrange.
   */
  typedef gMatrix<T_Value, SLICE> T_Rowsub;
  /**
   * Contiguous column subrange.
   */
  typedef const gMatrix<T_Value, COL_DENSE> T_ConstColsub;
  /**
   * Contiguous column subrange.
   */
  typedef gMatrix<T_Value, COL_DENSE> T_Colsub;
  /**
   * Const slice
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstSlice;
  /**
   * Slice
   */
  typedef gMatrix<T_Value, SLICE> T_Slice;
  /**
   * Dense type
   */
  typedef gMatrix<T_Value, COL_DENSE> T_Dense;
  /**
   * Const transpose
   */
  typedef const gMatrix<T_Value, ROW_DENSE> T_ConstTranspose;
  /**
   * Transpose
   */
  typedef gMatrix<T_Value, ROW_DENSE> T_Transpose;
}; // gMatrixTraits


/**
 * Matrix core class spcialization for dense col-major matrices (COL_DENSE)
 */
template<typename TF_Value>
class gMatrixCore< gMatrix<TF_Value, COL_DENSE> > : public gMatrixMemBase< gMatrix<TF_Value, COL_DENSE> > {
public:
  typedef gMatrix<TF_Value, COL_DENSE> T_Matrix;
  typedef gMatrixMemBase<T_Matrix> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_Storage T_Storage;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_ConstReference T_ConstReference;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_ConstSlice T_ConstSlice;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_ConstRow T_ConstRow;
  typedef typename T_Base::T_Row T_Row;
  typedef typename T_Base::T_ConstCol T_ConstCol;
  typedef typename T_Base::T_Col T_Col;
  typedef typename T_Base::T_ConstRowIterator T_ConstRowIterator;
  typedef typename T_Base::T_RowIterator T_RowIterator;
  typedef typename T_Base::T_ConstColIterator T_ConstColIterator;
  typedef typename T_Base::T_ColIterator T_ColIterator;
  typedef typename T_Base::T_ConstRowsub T_ConstRowsub;
  typedef typename T_Base::T_Rowsub T_Rowsub;
  typedef typename T_Base::T_ConstColsub T_ConstColsub;
  typedef typename T_Base::T_Colsub T_Colsub;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;
  typedef typename T_Base::T_StorageElem T_StorageElem;
  typedef typename T_Base::T_Dense T_Dense;

public:
  // no init ctor
  inline gMatrixCore(gLinearNoInit);

  // default ctor
  inline gMatrixCore(void);

  // ctor specifying size
  inline gMatrixCore(T_Index rows,
                     T_Index cols);

  // copy ctor
  inline gMatrixCore(const T_Matrix& mat);

  // ctor from dense gRowVector
  inline gMatrixCore(const gRowVector<T_Value, DENSE>& vec);

  // ctor from dense gColVector
  inline gMatrixCore(const gColVector<T_Value, DENSE>& vec);

  /**
   * Initializing ctor 2: use existing space read-only.
   */
  gMatrixCore(T_ConstPointer new_startAddr, T_Index rows, T_Index cols);

  // configure with supplied parameters
  inline void configure(T_Index new_nRows,
                        T_Index new_nCols,
                        const T_StorageElem* new_stor,
                        T_ConstPointer new_startAddr);

  // map const address from matrix coords
  inline T_Pointer mapAddr(T_Index row,
                           T_Index col) const;

  // dump info
  inline void specialInfo(int i) const;

  // transpose
  inline T_Transpose specialTranspose(void) const;

  // return row
  inline T_Row specialRow(T_Index row) const;

  // return col
  inline T_Col specialCol(T_Index col) const;

  // return row begin iterator
  inline T_RowIterator specialRowBegin(void) const;

  // return row end iterator
  inline T_RowIterator specialRowEnd(void) const;

  // return col begin iterator
  inline T_ColIterator specialColBegin(void) const;

  // return col end iterator
  inline T_ColIterator specialColEnd(void) const;

  // return reference to element
  inline T_Reference specialElem(T_Index row,
                                 T_Index col) const;

  // return dense row subrange
  inline T_Rowsub specialExtractRows(T_Index fromRow,
                                     T_Index toRow) const;

  // return dense col subrange
  inline T_Colsub specialExtractCols(T_Index fromCol,
                                     T_Index toCol) const;

  // return slice
  inline T_Slice specialSlice(const gIndexRange& row,
                              const gIndexRange& col) const;
}; // gMatrixCore


/**
 * Traits class for slice matrix (SLICE)
 */
template<typename TF_Value>
class gMatrixTraits< gMatrix<TF_Value, SLICE> > {
public:
  /**
   * Type of matrix.
   */
  typedef gMatrix<TF_Value, SLICE> T_Matrix;
  /**
   * Type of matrix elements.
   */
  typedef TF_Value T_Value;
  /**
   * Type of const (immuting) pointer to elements.
   */
  typedef const T_Value* T_ConstPointer;
  /**
   * Type of pointer to elements.
   */
  typedef T_Value* T_Pointer;
  /**
   * Type of const (immuting) reference to elements.
   */
  typedef const T_Value& T_ConstReference;
  /**
   * Type of reference to elements.
   */
  typedef T_Value& T_Reference;
  /**
   * Index type.
   */
  typedef typename gLinearTypes::T_Index T_Index;
  /**
   * Type of storage the implements.
   */
  typedef SLICE T_Storage;
  /**
   * Const row
   */
  typedef const gRowVector<T_Value, SLICE> T_ConstRow;
  /**
   * Row
   */
  typedef gRowVector<T_Value, SLICE> T_Row;
  /**
   * Const column
   */
  typedef const gColVector<T_Value, SLICE> T_ConstCol;
  /**
   * Column
   */
  typedef gColVector<T_Value, SLICE> T_Col;
  /**
   * Const row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_ConstRowIterator;
  /**
   * Row iterator
   */
  typedef gMatrixIterator<T_Matrix, ROW> T_RowIterator;
  /**
   * Const column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ConstColIterator;
  /**
   * Column iterator
   */
  typedef gMatrixIterator<T_Matrix, COL> T_ColIterator;
  /**
   * Const contiguous row subrange.
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstRowsub;
  /**
   * Contiguous row subrange.
   */
  typedef gMatrix<T_Value, SLICE> T_Rowsub;
  /**
   * Const contiguous column subrange.
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstColsub;
  /**
   * Contiguous column subrange.
   */
  typedef gMatrix<T_Value, SLICE> T_Colsub;
  /**
   * Const slice
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstSlice;
  /**
   * Slice
   */
  typedef gMatrix<T_Value, SLICE> T_Slice;
  /**
   * Dense type
   */
  typedef gMatrix<T_Value, ROW_DENSE> T_Dense;
  /**
   * Const transpose
   */
  typedef const gMatrix<T_Value, SLICE> T_ConstTranspose;
  /**
   * Transpose
   */
  typedef gMatrix<T_Value, SLICE> T_Transpose;
}; // gMatrixTraits

/**
 * Matrix core class spcialization for slice matrices (SLICE)
 */
template<typename TF_Value>
class gMatrixCore< gMatrix<TF_Value, SLICE> > : public gMatrixMemBase< gMatrix<TF_Value, SLICE> > {
public:
  typedef gMatrix<TF_Value, SLICE> T_Matrix;
  typedef gMatrixMemBase<T_Matrix> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_Storage T_Storage;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_ConstReference T_ConstReference;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_ConstSlice T_ConstSlice;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_ConstRow T_ConstRow;
  typedef typename T_Base::T_Row T_Row;
  typedef typename T_Base::T_ConstCol T_ConstCol;
  typedef typename T_Base::T_Col T_Col;
  typedef typename T_Base::T_ConstRowIterator T_ConstRowIterator;
  typedef typename T_Base::T_RowIterator T_RowIterator;
  typedef typename T_Base::T_ConstColIterator T_ConstColIterator;
  typedef typename T_Base::T_ColIterator T_ColIterator;
  typedef typename T_Base::T_ConstRowsub T_ConstRowsub;
  typedef typename T_Base::T_Rowsub T_Rowsub;
  typedef typename T_Base::T_ConstColsub T_ConstColsub;
  typedef typename T_Base::T_Colsub T_Colsub;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;
  typedef typename T_Base::T_StorageElem T_StorageElem;
  typedef typename T_Base::T_Dense T_Dense;

private:
  // row stride
  T_Index rowStride;
  // column stride
  T_Index colStride;

public:
  // no init ctor
  inline gMatrixCore(gLinearNoInit);

  // default ctor
  inline gMatrixCore(void);

  // ctor specifying size
  inline gMatrixCore(T_Index rows,
                     T_Index cols);

  // copy ctor
  inline gMatrixCore(const T_Matrix& mat);

  // init ctor
  inline gMatrixCore(T_Index new_nRows,
                     T_Index new_nCols,
                     const T_StorageElem* new_stor,
                     T_ConstPointer new_startAddr,
                     T_Index new_rowStride,
                     T_Index new_colStride);

  // ctor from dense gRowVector
  inline gMatrixCore(const gRowVector<T_Value, DENSE>& vec);

  // ctor from slice gRowVector
  inline gMatrixCore(const gRowVector<T_Value, SLICE>& vec);

  // ctor from dense gColVector
  inline gMatrixCore(const gColVector<T_Value, DENSE>& vec);

  // ctor from slice gColVector
  inline gMatrixCore(const gColVector<T_Value, SLICE>& vec);

  // configure with supplied parameters
  inline void configure(T_Index new_nRows,
                        T_Index new_nCols,
                        const T_StorageElem* new_stor,
                        T_ConstPointer new_startAddr,
                        T_Index new_rowStride,
                        T_Index new_colStride);

  // map const address from matrix coords
  inline T_Pointer mapAddr(T_Index row,
                           T_Index col) const;

  // dump info
  inline void specialInfo(int i) const;

  // transpose
  inline T_Transpose specialTranspose(void) const;

  // return row
  inline T_Row specialRow(T_Index row) const;

  // return col
  inline T_Col specialCol(T_Index col) const;

  // return row begin iterator
  inline T_RowIterator specialRowBegin(void) const;

  // return row end iterator
  inline T_RowIterator specialRowEnd(void) const;

  // return col begin iterator
  inline T_ColIterator specialColBegin(void) const;

  // return col end iterator
  inline T_ColIterator specialColEnd(void) const;

  // return reference to element
  inline T_Reference specialElem(T_Index row,
                                 T_Index col) const;

  // return dense row subrange
  inline T_Rowsub specialExtractRows(T_Index fromRow,
                                     T_Index toRow) const;

  // return dense col subrange
  inline T_Colsub specialExtractCols(T_Index fromCol,
                                     T_Index toCol) const;

  // return slice
  inline T_Slice specialSlice(const gIndexRange& row,
                              const gIndexRange& col) const;
}; // gMatrixCore


GLINEAR_NAMESPACE_END

#include "gMatrixCore.tcc"

#endif // GMATRIX_CORE_H
