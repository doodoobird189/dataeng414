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

#ifndef GVECTOR_MEM_BASE_H
#define GVECTOR_MEM_BASE_H

// glinear headers
#include "gLinearMem.h"
#include "gVectorBase.h"
#include "gLinearMacros.h"

GLINEAR_NAMESPACE_BEGIN

/** tag type to indicate the no-copy constructor */
struct vector_no_copy_tag { };

/**
 * Base class for all vectors storing data in memory.
 */
template<typename TF_Vector>
class gVectorMemBase : public gVectorBase<TF_Vector> {
public:
  /* types have to be fully qualified */
  typedef TF_Vector T_Vector;
  typedef gVectorBase<T_Vector> T_Base;
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

  typedef RefCntStorage<T_Value> T_StorageElem;
  friend class gVectorBase<T_Vector>;

  /**
   * Class name string
   */
  static const char className[];

private:
  /**
   * Number of elements.
   */
  T_Index nElem;
  /**
   * Pointer to storage class.
   */
  T_StorageElem* stor;
  /**
   * Address of the first element. Used to save extra indirection to storage
   * class.
   */
  T_Pointer startAddr;

public:
  /**
   * No init. Do nothing.
   */
  gVectorMemBase(gLinearNoInit);
  /**
   * Default constructor and constructor specifying size. \\
   * Constructs uninitialized gVectorMemBase containing {\em elem} elements.
   * A new storage element containing {\em elem} elements is allocated
   * and the reference count is set to 1. A default size of 0 is assumed
   * if this is called without parameters (as a default constructor).
   * @param elem Number of elements. Default = 0.
   */
  gVectorMemBase(T_Index new_nElem = 0);
  /**
   * Copy constructor. \\
   * Construct from another gVectorMemBase.
   * A {\bf shallow} copy is performed and the new gVectorMemBase is constructed
   * as an alias of the argument gVectorMemBase {\bf vec}. The reference count
   * of the storage class is incremented.
   * @param vec Source vector.
   */
  gVectorMemBase(const gVectorMemBase<T_Vector>& vec);
  /**
   * Construct with explicit initialization. \\
   * Construct new gVectorMemBase using the supplied storage class object,
   * start address and number of elements. Increment the reference
   * count on the storage class object.
   * @param new_nElem Number of elements (nElem).
   * @param new_stor RefCntStorage<T> storage class object (stor).
   * @param new_startAddr Address of first element (startAddr).
   */
  gVectorMemBase(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr);
  /**
   * Construct and initialize from existing C/C++ array. \\
   * A starting address {\em addr} and the number of elements (size of
   * vector), {\em elem} is specified. It is assumed that this memory
   * is already allocated and is readable.
   * New storage is allocated to contain {\em elem} T objects.
   * This new storage is then initialized through the assignment
   * operator from the {\em elem} elements stored at starting address
   * {\em addr}.
   * @param addr Starting address of contiguous block of elements.
   * @param elem Number of elements.
   */
  gVectorMemBase(T_ConstPointer addr,
                 T_Index elem);
  /**
   * Initializing ctor 2: use existing space read-only.
   */
  gVectorMemBase(T_ConstPointer addr,
                 T_Index elem,
                 const vector_no_copy_tag&);
  /*
   * Construct and initialize with supplied value. \\
   * The number of elements (size) and object with which to initialize
   * the vector is specified.
   * New storage capable of containing {\em elem} elements is allocated
   * and {\bf each} element is initialized to the argument {\em fill} through
   * the assignment operator.
   * @param fill Object to use as initializer.
   * @param elem Size (number of elements) of vector.
   */
  gVectorMemBase(T_ConstReference fill,
                 T_Index elem);
  /**
   * Destructor. \\
   * Decrement reference count of storage class and destruct it if
   * advisable. Advisability is determined from the return value of the
   * {\em T_StorageElem::rmref()} function.
   */
  ~gVectorMemBase(void);

  /*
   * Configure with explicit arguments. \\
   * Configure vector using the supplied storage class object,
   * start address and number of elements. Increment the reference
   * count of the supplied storage element.
   * @param new_nElem Number of elements (nElem).
   * @param new_stor RefCntStorage<T> storage class object (stor).
   * @param new_startAddr Address of first element (startAddr).
   */
  void configure(T_Index new_nElem,
                 const T_StorageElem* new_stor,
                 T_ConstPointer new_startAddr);

  /**
   * Dump info
   */
  void specialInfo(int i) const;
  /**
   * Return nElem
   */
  T_Index specialSize(void) const;

  /**
   * Set nElem
   */
  bool specialResize(T_Index new_nElem);

  /**
   * Return the number of references on this vector.
   */
  T_Index specialGetRefs(void) const;

  /**
   * Ensure that no aliases exist.
   */
  void specialDeepen(void) const;
  /**
   * Make this alias of other vector.
   */
  void specialAlias(const gVectorMemBase<T_Vector>& vec) const;
  /**
   * Make this an alias of other vector.
   */
  void specialAlias(gVectorMemBase<T_Vector>& vec);

  /*
   * Stuff for the world.
   */

  /**
   * Return the value of nElem.
   */
  T_Index getNElem(void) const;
  /**
   * Return a pointer to the T_StorageElem storage element.
   * @return Pointer to stor.
   */
  T_StorageElem* getStor(void) const;
  /**
   * Return a pointer to the first element of the vector (start address).
   * @return Pointer to the first element.
   */
  T_Pointer getStartAddr(void) const;

}; // gVectorMemBase

GLINEAR_NAMESPACE_END

#include "gVectorMemBase.tcc"

#endif // GVECTOR_MEM_BASE_H
