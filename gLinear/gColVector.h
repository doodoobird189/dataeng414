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

#ifndef GCOLVECTOR_H
#define GCOLVECTOR_H

// glinear headers
#include "gVectorCore.h"

GLINEAR_NAMESPACE_BEGIN

// forward decl.
template<typename T_Value, typename T_Storage>
class gRowVector;

template<typename T_Value, typename T_Storage = DENSE>
class gColVector;

/**
 * Traits class specializations for a vector using DENSE storage.
 */
template<typename TF_Value>
class gVectorTraits< gColVector<TF_Value, DENSE> > {
public:
  /**
   * Type of vector.
   */
  typedef gColVector<TF_Value, DENSE> T_Vector;
  /**
   * Type of vector elements.
   */
  typedef TF_Value T_Value;
  /**
   * Type of const (immuting) pointer to vector elements.
   */
  typedef const T_Value* T_ConstPointer;
  /**
   * Type of pointer to vector elements.
   */
  typedef T_Value* T_Pointer;
  /**
   * Type of const (immuting) reference to vector elements.
   */
  typedef const T_Value& T_ConstReference;
  /**
   * Type of reference to vector elements.
   */
  typedef T_Value& T_Reference;
  /**
   * Index type.
   */
  typedef typename gLinearTypes::T_Index T_Index;
  /**
   * Type of storage the vector implements.
   */
  typedef DENSE T_Storage;
  /**
   * Type of const iterator used to traverse the vector.
   */
  typedef const T_Value* T_ConstIterator;
  /**
   * Type of iterator used to traverse the vector.
   */
  typedef T_Value* T_Iterator;
  /**
   * Const subrange of vector.
   */
  typedef const gColVector<T_Value, DENSE> T_ConstSubrange;
  /**
   * Subrange of vector.
   */
  typedef gColVector<T_Value, DENSE> T_Subrange;
  /**
   * Const slice of vector.
   */
  typedef const gColVector<T_Value, SLICE> T_ConstSlice;
  /**
   * Slice of vector.
   */
  typedef gColVector<T_Value, SLICE> T_Slice;
  /**
   * Dense of vector.
   */
  typedef gColVector<T_Value, DENSE> T_Dense;
  /**
   * Const transpose of vector.
   */
  typedef const gRowVector<T_Value, DENSE> T_ConstTranspose;
  /**
   * Transpose of vector.
   */
  typedef gRowVector<T_Value, DENSE> T_Transpose;
}; // gVectorTraits

/**
 * Column vector implementing DENSE storage.
 */
template<typename TF_Value>
class gColVector<TF_Value, DENSE> : public gDenseVector< gColVector<TF_Value, DENSE> > {
public:
  typedef gColVector<TF_Value, DENSE> T_Vector;
  typedef gDenseVector<T_Vector> T_Base;
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
  typedef typename T_Base::T_StorageElem T_StorageElem;
  typedef typename T_Base::T_Dense T_Dense;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;
  /**
   * Class name string
   */
  static const char className[];

  /**
   * No init ctor.
   */
  gColVector(gLinearNoInit);

  /**
   * Default ctor and ctor specifying size.
   */
  explicit gColVector(T_Index elem = 0);

  /**
   * Copy ctor.
   */
  gColVector(const gColVector<TF_Value, DENSE>& vec);

  /**
   * Initializing ctor.
   */
  gColVector(T_Index new_nElem,
             const T_StorageElem* new_stor,
             T_ConstPointer new_startAddr);

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
  gColVector(T_ConstPointer new_startAddr, T_Index len);

  /**
   * Copy assignment.
   *
   * Needed because the operator= doesn't inherit and we don't want
   * the compiler to create an implicit one.
   */
  T_Vector& operator=(const T_Vector& vec);

  /** gColVector assignment. Needed because the operator= doesn't inherit.
   */
  template<typename T_SrcValue, typename T_SrcStorage>
  T_Vector& operator=(const gColVector<T_SrcValue, T_SrcStorage>& src);

  void specialInfo(int i) const;

}; // gColVector

/**
 * Traits class specializations for a vector using SLICE storage.
 */
template<typename TF_Value>
class gVectorTraits< gColVector<TF_Value, SLICE> > {
public:
  /**
   * Type of vector.
   */
  typedef gColVector<TF_Value, SLICE> T_Vector;
  /**
   * Type of vector elements.
   */
  typedef TF_Value T_Value;
  /**
   * Type of const (immuting) pointer to vector elements.
   */
  typedef const T_Value* T_ConstPointer;
  /**
   * Type of pointer to vector elements.
   */
  typedef T_Value* T_Pointer;
  /**
   * Type of const (immuting) reference to vector elements.
   */
  typedef const T_Value& T_ConstReference;
  /**
   * Type of reference to vector elements.
   */
  typedef T_Value& T_Reference;
  /**
   * Index type.
   */
  typedef typename gLinearTypes::T_Index T_Index;
  /**
   * Type of storage the vector implements.
   */
  typedef SLICE T_Storage;
  /**
   * Type of const iterator used to traverse the vector.
   */
  typedef gSliceVectorIterator<T_Value> T_ConstIterator;
  /**
   * Type of iterator used to traverse the vector.
   */
  typedef gSliceVectorIterator<T_Value> T_Iterator;
  /**
   * Const subrange of vector.
   */
  typedef const gColVector<T_Value, SLICE> T_ConstSubrange;
  /**
   * Subrange of vector.
   */
  typedef gColVector<T_Value, SLICE> T_Subrange;
  /**
   * Const slice of vector.
   */
  typedef const gColVector<T_Value, SLICE> T_ConstSlice;
  /**
   * Slice of vector.
   */
  typedef gColVector<T_Value, SLICE> T_Slice;
  /**
   * Dense of vector.
   */
  typedef gColVector<T_Value, DENSE> T_Dense;
  /**
   * Const transpose of vector.
   */
  typedef const gRowVector<T_Value, DENSE> T_ConstTranspose;
  /**
   * Transpose of vector.
   */
  typedef gRowVector<T_Value, SLICE> T_Transpose;
}; // gVectorTraits

/**
 * Column vector implementing SLICE storage.
 */
template<typename TF_Value>
class gColVector<TF_Value, SLICE> : public gSliceVector< gColVector<TF_Value, SLICE> > {
public:
  typedef gColVector<TF_Value, SLICE> T_Vector;
  typedef gSliceVector<T_Vector> T_Base;
  typedef typename T_Base::T_Value T_Value;
  typedef typename T_Base::T_ConstPointer T_ConstPointer;
  typedef typename T_Base::T_Pointer T_Pointer;
  typedef typename T_Base::T_ConstReference T_ConstReference;
  typedef typename T_Base::T_Reference T_Reference;
  typedef typename T_Base::T_Slice T_Slice;
  typedef typename T_Base::T_Subrange T_Subrange;
  typedef typename T_Base::T_Index T_Index;
  typedef typename T_Base::T_ConstIterator T_ConstIterator;
  typedef typename T_Base::T_Iterator T_Iterator;
  typedef typename T_Base::T_StorageElem T_StorageElem;
  typedef typename T_Base::T_Dense T_Dense;
  typedef typename T_Base::T_ConstTranspose T_ConstTranspose;
  typedef typename T_Base::T_Transpose T_Transpose;

  /**
   * Class name string
   */
  static const char className[];

  /**
   * No init ctor.
   */
  gColVector(gLinearNoInit);

  /**
   * Default ctor and ctor specifying size.
   */
  explicit gColVector(T_Index elem = 0);

  /**
   * Copy ctor.
   */
  gColVector(const gColVector<TF_Value, SLICE>& vec);

  /**
   * Initializing ctor.
   */
  gColVector(T_Index new_nElem,
             const T_StorageElem* new_stor,
             T_ConstPointer new_startAddr,
             T_Index new_stride);

  /** Copy assignment. Needed because the operator= doesn't inherit and we don't want the
   * compiler to create an implicit one.
   */
  T_Vector& operator=(const T_Vector& vec);

  /** gColVector assignment. Needed because the operator= doesn't inherit.
   */
  template<typename T_SrcValue, typename T_SrcStorage>
  T_Vector& operator=(const gColVector<T_SrcValue, T_SrcStorage>& src);

  void specialInfo(int i) const;

}; // gColVector

GLINEAR_NAMESPACE_END

#include "gColVector.tcc"

#endif // GCOLVECTOR_H
