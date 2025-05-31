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

#ifndef GLINEAR_MEM_H
#define GLINEAR_MEM_H

// glinear headers
#include "gError.h"
#include "gLinearMacros.h"
#include "gLinearAlloc.h"

GLINEAR_NAMESPACE_BEGIN

/**
 * Low-level one-dimensional non-growable homogeneous container class
 * that implements simple reference counting.
 * Increasing or decreasing the reference count
 * is done explicitly (through the addref() and rmref() functions).
 * This class manages a homogenous container of fixed size.
 */
template<typename T>
class RefCntStorage {
private:
  /**
   * Type of elements that this stores.
   */
  typedef T T_Value;
  /**
   * Pointer to element type.
   */
  typedef T_Value* T_Pointer;
  /**
   * Reference to element type.
   */
  typedef T_Value& T_Reference;
  /**
   * Index type
   */
  typedef gLinearTypes::T_Index T_Index;

  static const char ClassName[];
  /**
   * Number of class instances referring to this storage class. \\
   * Should normally be integral positive non-zero. (May be zero after
   * invocation of default ctor).
   * A value of zero indicates that it would do no harm to destruct
   * this class (it is in fact advisable in order to prevent memory leaks).
   */
  T_Index refs;
  /**
   * Size of the container. \\
   * Number of elements (of type T) that this container may contain in
   * the allocated storage.
   * Should be integral zero or positive.
   */
  T_Index size;
  /**
   * Starting address of storage. \\
   * Storage is allocated as a C++ array ( T[] ). It is up to the class
   * referencing this storage to attach any higher-level interpretation
   * to specific locations in this array. \\
   * <bitch> \\
   * Just to ensure that I make myself clear: This class does not know
   * or (more importantly) care about rows, columns, transposes,
   * dense, sparse or for that matter matrices or vectors. It's only purpose
   * in life is to encapsulate the underlying storage and provide us with
   * the reference counting.
   * </bitch>
   */
  T_Pointer startAddr;
public:
  /**
   * Default contructor. \\
   * Construct with some sane (mostly zero) default values. \\
   * \begin{itemize}
   * \item No storage is allocated.
   * \item The reference count is set equal to zero.
   * \item The size is set equal to zero.
   * \end{itemize}
   */
  RefCntStorage(void);

  /**
   * Constructor specifying number of elements. \\
   * \begin{itemize}
   * \item An array of {\em newsize} T objects is allocated.
   * \item The reference count is set equal to zero.
   * \item The size is set to {\em newsize}.
   * \end{itemize}
   */
  RefCntStorage(T_Index new_size);

  /**
   * Construct and initialize from C/C++ array. \\
   * \begin{itemize}
   * \item Allocate array of {\em newsize} elements.
   * \item Set reference count equal to zero.
   * \item Set size to {\em newsize}.
   * \item Copy {\em newsize} elements from arrray starting at {\em addr}
   * into allocated storage.
   * \end{itemize}
   */
  RefCntStorage(const T_Pointer addr, T_Index new_size);

  /**
   * Construct and initialize with supplied argument. \\
   * \begin{itemize}
   * \item Allocate array of {\em newsize} elements.
   * \item Set reference count equal to zero.
   * \item Set size to {\em newsize}.
   * \item Set each element of the allocated storage equal to {\em fill}.
   * \end{itemize}
   */
  RefCntStorage(const T_Value& fill, T_Index new_size);

  /**
   * Destructor. \\
   * Deallocate storage (delete array) if any is allocated.
   */
  ~RefCntStorage(void);

  /**
   * Dump info.
   */
  void info(int i) const;

  /**
   * Decrement reference count. \\
   * Decrements the reference count by one and deallocate {\em this} object
   * if there are no more references to it.
   * A {\em true} return value indicates that the {\em this}
   * object is unaware of any more references to it and {\em already} destructed itself
   * (in a good way that is).
   * @return {\bf true -} No more objects referencing this storage. Please
   * do not refererence it. \\
   * {\bf false -} One or more objects are still referencing this storage.
   * DONT DESTRUCT.
   */
  bool rmref(void);

  /**
   * Increment reference count.
   */
  void addref(void);

  /**
   * Return the number of storage elements that this object contains
   * (allocated).\\
   * @return Number of storage elements. Positive or zero integer.
   */
  T_Index getSize(void);

  /**
   * Change the size (number of allocated storage elements) that this
   * object contains. \\
   * \begin{itemize}
   * \item The current storage (array) is first deallocated and
   * new storage is then reallocated.
   * \item The {\em size} is changed to {\em newsize}.
   * \item The reference count is left untouched.
   * \end{itemize}
   * This basically has the same effect as calling
   * RefCntStorage(T_Index newsize), but retaining the reference counting
   * values. \\
   * <bitch> \\
   * Please save yourself a lot of misery and do not make any assumptions
   * about your data still being intact after calling this function. \\
   * </bitch>
   */
  void setSize(T_Index new_size);

  /**
   * Return the number of references to this object.
   * @return Number of references to this object.
   */
  T_Index getRefs(void);

  /**
   * Return address to start of storage (array). \\
   * <bitch> \\
   * {\bf Note: } This function is only included to allow for some
   * optimization in the matrix and vector classes using this class. \\
   * If you save the address returned by this function and
   * dereference it {\em after} destructing this class you are {\bf on your
   * own}. Anything might happen. You have been warned. Have fun.\\
   * </bitch>
   */
  T* getStartAddr(void);

  /**
   * Make an independent duplicate of this object. \\
   * A new RefCntStorage<T> object, of the same size as this one, is
   * created.
   * The elements are then initialized by performing a
   * deep copy of the data (pointed to by start_addr) from this object
   * to the new one. The reference count is set to 0.
   * A pointer to the new RefCntStorage<T> is returned.
   * This is the only chance you have of saving this pointer. It is not
   * automatically saved {\bf anywhere}.\\
   * Call this a few million times wihtout saving the returned pointer
   * if you are a fan of huge memory leaks.
   * @return Pointer to the RefCntStorage<T> object that is a duplicate
   * of this one.
   */
  RefCntStorage<T_Value>* duplicate(void);

  /**
   * Set the RefCntStorage to an existing array of the right
   * type. This method is very special: it makes sets the refcount to
   * 2.  To destroy this object, you *must* call release() beforehand,
   * otherwise it will attempt to delete[] its storage.
   *
   * WARNING: Make very sure that you only pass this on to const
   * objects.  It may not ever be modified.
   */
  void setExisting(const T* stor, int len);

  /**
   * invalidates this object, and makes sure that it does not delete
   * its storage.
   */
  void releaseExisting();

  /**
   * Custom new.
   */
  void* operator new(size_t size);

  /**
   * Custom delete.
   */
  void operator delete(void* stor);

}; // RefCntStorage

GLINEAR_NAMESPACE_END

#include "gLinearMem.tcc"

#endif // GLINEAR_MEM_H
