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

#ifndef GMATRIX_H
#define GMATRIX_H

// glinear headers
#include "gMatrixCore.h"

GLINEAR_NAMESPACE_BEGIN

/**
 * Matrix class. \\
 * This contains all the members that doesn't normally inherit,
 * i.e. the constuctors and the assignment operators.
 */
template<typename TF_Value, typename TF_Storage = ROW_DENSE>
class gMatrix : public gMatrixCore< gMatrix<TF_Value, TF_Storage> > {
public:
  typedef gMatrix<TF_Value, TF_Storage> T_Matrix;
  typedef gMatrixCore<T_Matrix> T_Base;
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

  /**
   * No Init constuctor. \\
   */
  gMatrix(gLinearNoInit);

  /**
   * Default constructor. \\
   * Construct [0 x 0] (empty) matrix.
   */
  gMatrix(void);

  /**
   * Constructor specifying size. \\
   * Construct [{\em rows x cols}] uninitialised matrix.
   * @param rows Number of rows.
   * @param cols Number of cols.
   */
  gMatrix(T_Index rows,
          T_Index cols);

  /**
   * Copy contructor. \\
   * Construct as a shallow copy from the argument matrix.
   * @param mat Argument matrix.
   */
  gMatrix(const T_Matrix& mat);

  /**
   * Construct from gRowVector
   */
  template<typename T_VecStor>
  explicit gMatrix(const gRowVector<T_Value, T_VecStor>& vec);

  /**
   * Construct from gColVector
   */
  template<typename T_VecStor>
  explicit gMatrix(const gColVector<T_Value, T_VecStor>& vec);

  /**
   * Initializing ctor 2: Careful with this one. Of course the memory
   * new_startAddr points to must really be allocated/available. But
   * also, lifetime control of this memory is entirely up to whoever
   * allocated it originally. If for some reason it disappears before
   * this object does, well, then you are in trouble.
   *
   * ODD: I dont quite get the ConstPointer bit, it seems that we can
   * (and should) read and write to an object constructed in this
   * manner.
   *
   * TODO: We should at some point write a proper implementation of
   * this that takes a proper object as the first parameter so that we
   * can make use of refcounting.
   */
  gMatrix(T_ConstPointer new_startAddr, T_Index rows, T_Index cols);

  /**
   * Copy assignment. \\
   * Doesn't inherit.
   */
  gMatrix& operator=(const gMatrix& mat);

  /**
   * gMatrix assignment
   */
  template<typename T_SrcValue, typename T_SrcStorage>
  gMatrix& operator=(const gMatrix<T_SrcValue, T_SrcStorage>& mat);

  void specialInfo(int i) const;

}; // gMatrix

template<typename T>
using gRowMatrix = gMatrix<T, ROW_DENSE>;

template<typename T>
using gColMatrix = gMatrix<T, COL_DENSE>;


GLINEAR_NAMESPACE_END

#include "gMatrix.tcc"

#endif // GMATRIX_H
