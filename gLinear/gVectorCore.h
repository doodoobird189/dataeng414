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

#ifndef GVECTOR_CORE_H
#define GVECTOR_CORE_H

// glinear headers
#include "gVectorIterator.h"
#include "gVectorMemBase.h"
#include "gLinearMacros.h"

GLINEAR_NAMESPACE_BEGIN

/*
 * Generic Vector class that is available in a number of implementations.
 * This class is specialized on T_Storage to create different types (e.g. slice,
 * dense, sparse) of vectors.
 */
template<typename T_Vector>
class gDenseVector;

template<typename T_Vector>
class gSliceVector;

/**
 * Generic vector implementing dense storage.
 */
template<typename TF_Vector>
class gDenseVector : public gVectorMemBase<TF_Vector> {
public:
  typedef TF_Vector T_Vector;
  typedef gVectorMemBase<T_Vector> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_ConstReference T_ConstReference;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_ConstSlice T_ConstSlice;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_ConstSubrange T_ConstSubrange;
  typedef typename T_Base::T_Subrange T_Subrange;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstIterator T_ConstIterator;
  typedef typename T_Base::T_Iterator T_Iterator;
  typedef typename T_Base::T_Dense T_Dense;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;

  typedef typename gVectorMemBase<T_Vector>::T_StorageElem T_StorageElem;
  friend class gVectorBase<T_Vector>;

  // class name string
  static const char className[];

  // no init ctor.
  gDenseVector(gLinearNoInit);

  // default ctor. and ctor. specifying size.
  gDenseVector(T_Index elem = 0);

  // copy ctor.
  gDenseVector(const gDenseVector<T_Vector>& vec);

  // initializing ctor.
  gDenseVector(T_Index new_nElem,
               const T_StorageElem* new_stor,
               T_ConstPointer new_startAddr);

  /**
   * Initializing ctor 2: use existing space read-only.
   */
  gDenseVector(T_ConstPointer new_startAddr,
         T_Index len);

  // configure
  void configure(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr);

protected:
  /* Specialized implementations of functions
   */

  // dump info
  void specialInfo(int i) const;

  // element access
  T_Reference specialElem(T_Index elem) const;

  // begin iterator
  T_Iterator specialBegin(void) const;

  // past-the-end iterator
  T_Iterator specialEnd(void) const;

  // slice
  T_Slice specialSlice(const gIndexRange& elems) const;

  // subrange
  T_Subrange specialSubrange(T_Index fromElem,
                             T_Index toElem) const;

  // transpose
  T_Transpose specialTranspose(void) const;

}; // gDenseVector


/**
 * Generic vector implementing slice storage
 */
template<typename TF_Vector>
class gSliceVector : public gVectorMemBase<TF_Vector> {
public:
  typedef TF_Vector T_Vector;
  typedef gVectorMemBase<T_Vector> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_ConstSlice T_ConstSlice;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_ConstSubrange T_ConstSubrange;
  typedef typename T_Base::T_Subrange T_Subrange;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstIterator T_ConstIterator;
  typedef typename T_Base::T_Iterator T_Iterator;
  typedef typename T_Base::T_Dense T_Dense;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;

  typedef typename gVectorMemBase<T_Vector>::T_StorageElem T_StorageElem;
  friend class gVectorBase<T_Vector>;

  // class name string
  static const char className[];

private:
  /**
   * Stride of the vector.
   */
  T_Index stride;
protected:
  // no init ctor.
  gSliceVector(gLinearNoInit);

protected:
  // default ctor. and ctor. specifying size
  gSliceVector(T_Index elem = 0);

public:
  // copy ctor.
  gSliceVector(const gSliceVector<T_Vector>& vec);

  // init ctor.
  gSliceVector(T_Index new_nElem,
               const T_StorageElem* new_stor,
               T_ConstPointer new_startAddr,
               T_Index new_stride);

  // configure with supplied parameters
  void configure(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr,
                 T_Index new_stride);

  /* Required Specialized Implementations
   */
protected:
  // dump info
  void specialInfo(int i) const;

  // resize: illegal on slice unless resize does nothing (which
  // gets picked up in gVectorBase<T_Vector>
  bool specialResize(T_Index elem);

  // element access
  T_Reference specialElem(T_Index elem) const;

  // begin iterator
  T_Iterator specialBegin(void) const;

  // past-the-end iterator
  T_Iterator specialEnd(void) const;

  // slice
  T_Slice specialSlice(const gIndexRange& elems) const;

  // subrange
  T_Subrange specialSubrange(T_Index fromElem,
                             T_Index toElem) const;

  // transpose
  T_Transpose specialTranspose(void) const;

public:
  /** Return stride
   */
  T_Index getStride(void) const;

}; // gSliceVector

GLINEAR_NAMESPACE_END

#include "gVectorCore.tcc"

#endif // GVECTOR_CORE_H
