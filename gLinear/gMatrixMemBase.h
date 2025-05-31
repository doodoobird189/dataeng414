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

#ifndef GMATRIX_MEM_BASE_H
#define GMATRIX_MEM_BASE_H

// glinear headers
#include "gLinearMacros.h"
#include "gLinearMem.h"
#include "gMatrixBase.h"


GLINEAR_NAMESPACE_BEGIN

/** tag type to indicate the no-copy constructor */
struct matrix_no_copy_tag { };

/**
 * Base class for all matrices storing data in memory.
 */
template<typename TF_Matrix>
class gMatrixMemBase : public gMatrixBase<TF_Matrix> {
public:
  // typedefs brought from T_Base
  typedef TF_Matrix T_Matrix;
  typedef gMatrixBase<T_Matrix> T_Base;

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
  typedef typename T_Base::T_ConstRowsub T_ConstRowsub;
  typedef typename T_Base::T_Rowsub T_Rowsub;
  typedef typename T_Base::T_ConstColsub T_ConstColsub;
  typedef typename T_Base::T_Colsub T_Colsub;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;
  typedef typename T_Base::T_Dense T_Dense;

  // type of storage element
  typedef RefCntStorage<T_Value> T_StorageElem;

private:
  /**
   * Number of rows. \\
   * Zero or positive integer.
   */
  T_Index nRows;
  /**
   * Number of columns. \\
   * Zero or positive integer.
   */
  T_Index nCols;
  /**
   * Storage.
   */
  T_StorageElem* stor;
  /**
   * Address of first element;
   */
  T_Pointer startAddr;
protected:
  /**
   * No init. Do nothing.
   */
  gMatrixMemBase(gLinearNoInit);
  /**
   * Default constructor. \\
   * Consructs a [0 x 0] matrix.
   */
  gMatrixMemBase(void);
  /**
   * Constructor specifying size. \\
   * Constructs uninitialized [new_nRows x new_nCols] gMatrixMemBase.
   * A new storage element containing {\em new_nRows * new_nCols}
   * elements is allocated and the reference count is set to 1.
   * @param new_nRows Number of rows.
   * @param new_nCols Number of cols.
   */
  gMatrixMemBase(T_Index new_nRows,
                 T_Index new_nCols);
  /**
   * Copy constructor. \\
   * Construct from another gMatrixMemBase.
   * A {\bf shallow} copy is performed and the new gMatrixMemBase is
   * constructed as an alias of the argument gMatrixMemBase {\bf mat}.
   * The reference count of the storage class is incremented.
   * @param mat Source matrix.
   */
  gMatrixMemBase(const gMatrixMemBase<T_Matrix>& mat);
  /**
   * Construct with explicit initialization. \\
   * Construct new gMatrixMemBase using the supplied storage class object,
   * start address and number of elements. Increment the reference
   * count on the storage class object.
   * @param new_nElem Number of elements (nElem).
   * @param new_stor RefCntStorage<T> storage class object (stor).
   * @param new_startAddr Address of first element (startAddr).
   */
  gMatrixMemBase(T_Index new_nRows,
                 T_Index new_nCols,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr);

  /**
   * Construct and initialize from existing C/C++ array.
   *
   * @param addr Starting address of contiguous block of elements.
   * @param rows Number of rows.
   * @param cols Number of cols.
   */
  gMatrixMemBase(T_ConstPointer addr, T_Index rows, T_Index cols);

  /**
   * Initializing ctor 2: use existing space read-only.
   */
  gMatrixMemBase(T_ConstPointer addr,
                 T_Index rows, T_Index cols,
                 const matrix_no_copy_tag&);

  /**
   * Destructor. \\
   * Decrement reference count of storage class and destruct it if
   * advisable. Advisability is determined from the return value of the
   * {\em T_StorageElem::rmref()} function.
   */
  ~gMatrixMemBase(void);
  /**
   * Configure with explicit arguments. \\
   * Configure matrix using the supplied storage class object,
   * start address and number of rows and columns. Increment the reference
   * count of the supplied storage element.
   * @param new_nRows Number of rows (nRows).
   * @param new_nCols Number of columns (nCols).
   * @param new_stor T_StorageElem class pointer (stor).
   * @param new_startAddr Address of first element (startAddr).
   */
  void configure(T_Index new_nRows,
                 T_Index new_nCols,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr);

public:
  /**
   * Dump info
   */
  void specialInfo(int i) const;
  /**
   * Return nRows
   */
  T_Index specialRows(void) const;
  /**
   * Return nCols
   */
  T_Index specialCols(void) const;
  /**
   * Set nRows and nCols
   */
  bool specialResize(T_Index new_nRows,
                     T_Index new_nCols);
  /**
   * Return the number of references on this matrix.
   */
  T_Index specialGetRefs(void) const;
  /**
   * Ensure that no aliases exist.
   */
  void specialDeepen(void) const;
  /**
   * Make this alias of other matrix.
   */
  void specialAlias(const gMatrixMemBase<T_Matrix>& vec) const;
  /**
   * Make this an alias of other matrix.
   */
  void specialAlias(gMatrixMemBase<T_Matrix>& vec);

  /*
   * Support
   */
public:
  /**
   * Return the value of nRows
   */
  T_Index getNRows(void) const;
  /**
   * Return the value of nCols
   */
  T_Index getNCols(void) const;
  /**
   * Return a pointer to the T_StorageElem storage element.
   * @return Pointer to stor.
   */
  const T_StorageElem* getStor(void) const;
  /**
   * Return a const pointer to the T_StorageElem storage element.
   * @return Const pointer to stor.
   */
  T_StorageElem* getStor(void);
  /**
   * Return const pointer to the first element of the allocated
   * storage (start address).
   * @return Const pointer to the allocated storage.
   */
  T_ConstPointer getStartAddr(void) const;
  /**
   * Return pointer to the first element of the allocated storage
   * (start address).
   * @return Pointer to the allocated storage.
   */
  T_Pointer getStartAddr(void);
}; // gMatrixMemBase

GLINEAR_NAMESPACE_END

#include "gMatrixMemBase.tcc"

#endif // GMATRIX_MEM_BASE_H
