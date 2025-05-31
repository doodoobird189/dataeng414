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

#ifndef GVECTOR_ITERATOR_H
#define GVECTOR_ITERATOR_H

#include <iterator>

// glinear headers
#include "gLinearMacros.h"

// standard headers
#include <iterator>

GLINEAR_NAMESPACE_BEGIN

/**
 * Iterator used to traverse slice vectors and matrices.
 */
template<typename TF_Value>
class gSliceVectorIterator :
  public std::iterator<std::bidirectional_iterator_tag, TF_Value, gLinearTypes::T_Index> {
public:
  typedef typename gLinearTypes::T_Index T_Index;
  typedef TF_Value T_Value;
  typedef const T_Value* T_ConstPointer;
  typedef T_Value* T_Pointer;
  typedef const T_Value& T_ConstReference;
  typedef T_Value& T_Reference;

private:
  /**
   * Pointer to current element.
   */
  T_Pointer addr;
  /**
   * Stride.
   */
  T_Index stride;
public:
  /**
   * Default contructor.
   * Zero-initialize the iterator. Don't use iterator before a value has been assigned.
   */
  gSliceVectorIterator(void);
  /**
   * Construct and initialize.
   * @param new_addr Address of initial element.
   * @param new_stride The number of elements between each iteration.
   */
  gSliceVectorIterator(T_ConstPointer new_addr,
                       T_Index new_stride);
  /**
   * Copy constuctor.
   * Construct and return duplicate of argument.
   */
  gSliceVectorIterator(const gSliceVectorIterator<T_Value>& iter);
  /**
   * Assigment operator.
   */
  gSliceVectorIterator operator=(const gSliceVectorIterator<T_Value>& iter);
  /**
   * const Dereference.
   */
  T_ConstReference operator*(void) const;
  /**
   * Dereference.
   */
  T_Reference operator*(void);
  /**
   * Equality operator.
   */
  bool operator==(const gSliceVectorIterator<T_Value>& iter) const;
  /**
   * Inequality operator.
   */
  bool operator!=(const gSliceVectorIterator<T_Value>& iter) const;
  /**
   * Postfix increment.
   */
  gSliceVectorIterator<T_Value> operator++(int i);
  /**
   * Prefix increment.
   */
  gSliceVectorIterator<T_Value> operator++(void);
  /**
   * Postfix decrement.
   */
  gSliceVectorIterator<T_Value> operator--(int i);
  /**
   * Prefix decrement.
   */
  gSliceVectorIterator<T_Value> operator--(void);
  /**
   * const Addition.
   */
  const gSliceVectorIterator<T_Value> operator+(T_Index delta) const;
  /**
   * Addition.
   */
  gSliceVectorIterator<T_Value> operator+(T_Index delta);
  /**
   * const Subtraction.
   */
  const gSliceVectorIterator<T_Value> operator-(T_Index delta) const;
  /**
   * Subtraction.
   */
  gSliceVectorIterator<T_Value> operator-(T_Index delta);
}; // gSliceVectorIterator

GLINEAR_NAMESPACE_END

#include "gVectorIterator.tcc"

#endif // GVECTOR_ITERATOR_H
