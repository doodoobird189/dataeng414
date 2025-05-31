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

#ifndef GMATRIX_ITERATOR_H
#define GMATRIX_ITERATOR_H

GLINEAR_NAMESPACE_BEGIN

// forward decl.
template<typename T_Matrix, typename T_VectorTag>
class gMatrixIterator;

// primary decl.
template<typename T_Matrix, typename T_VectorTag>
class gMatrixIteratorSupport;

// row iterator spec. for support
template<typename T_Matrix>
class gMatrixIteratorSupport<T_Matrix, ROW> {
public:
  typedef gMatrixIterator<T_Matrix, ROW> T_Iter;
  typedef typename T_Matrix::T_Row T_Row;
  typedef T_Row T_Vector;
public:
  static T_Row deref(const T_Iter& iter);
}; // gMatrixIteratorSupport

// column iterator spec. for support
template<typename T_Matrix>
class gMatrixIteratorSupport<T_Matrix, COL> {
public:
  typedef gMatrixIterator<T_Matrix, COL> T_Iter;
  typedef typename T_Matrix::T_Col T_Col;
  typedef T_Col T_Vector;
public:
  static T_Col deref(const T_Iter& iter);
}; // gMatrixIteratorSupport

/**
 * A general matrix iterator that dereferences to a row or a column. \\
 * T_Matrix indicates the type of the matrix returns this iterator. \\
 * T_Vector is the type of vector returned when the iterator is dereferenced.
 */
template<typename T_Matrix, typename T_VectorTag>
class gMatrixIterator
{
public:
  typedef typename T_Matrix::T_Index T_Index;
  typedef gMatrixIterator<T_Matrix, T_VectorTag> T_This;
  typedef gMatrixIteratorSupport<T_Matrix, T_VectorTag> T_Support;
  typedef typename T_Support::T_Vector T_Vector;
  // allow the support funcs. to dig in
  friend class gMatrixIteratorSupport<T_Matrix, T_VectorTag>;

  // standard iterator traits: 24.3.3
  typedef T_Vector value_type;
  typedef T_Index difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef std::bidirectional_iterator_tag iterator_category;

private:
  // index
  T_Index indx;
  // matrix
  T_Matrix* matrix;

public:
  /**
   * Construct uninitialized iterator.
   */
  gMatrixIterator(void);
  /**
   * Construct and initialize iterator.
   */
  gMatrixIterator(T_Index new_indx,
                  const T_Matrix* new_matrix);
  /**
   * Copy constructor.
   */
  gMatrixIterator(const T_This& iter);

  /**
   * Copy assignment.
   */
  T_This& operator=(T_This& iter);

public:
  /**
   * Const dereference.
   */
  inline const T_Vector operator*(void) const;
  /**
   * Dereference.
   */
  inline T_Vector operator*(void);
  /**
   * Equality operator.
   */
  bool operator==(const T_This& iter);
  /**
   * Inequality operator.
   */
  bool operator!=(const T_This& iter);

protected:
  // postfix inc support
  T_This postfixInc(void);
  // prefix inc support
  T_This& prefixInc(void);
  // postfix dec support
  T_This postfixDec(void);
  // prefix dec support
  T_This& prefixDec(void);

public:
#if 0
  /**
   * Const postfix increment.
   */
  const T_This operator++(int) const;
#endif
  /**
   * Postfix increment.
   */
  T_This operator++(int);
#if 0
  /**
   * Const prefix increment.
   */
  const T_This& operator++(void) const;
#endif
  /**
   * Prefix increment.
   */
  T_This& operator++(void);
#if 0
  /**
   * Const postfix decrement.
   */
  const T_This operator--(int) const;
#endif
  /**
   * Postfix decrement.
   */
  T_This operator--(int);
#if 0
  /**
   * Const prefix decrement.
   */
  const T_This& operator--(void) const;
#endif
  /**
   * Prefix decrement.
   */
  T_This& operator--(void);
protected:
  // addition support
  T_This add(T_Index i) const;
  // subtraction support
  T_This subtract(T_Index i) const;
 public:
  /**
   * Const addition.
   */
  const T_This operator+(T_Index i) const;
  /**
   * Adition.
   */
  T_This operator+(T_Index i);
  /**
   * Const subtraction.
   */
  const T_This operator-(T_Index i) const;
  /**
   * Subtraction.
   */
  T_This operator-(T_Index i);

}; // gMatrixIterator

GLINEAR_NAMESPACE_END

#include "gMatrixIterator.tcc"

#endif // GMATRIX_ITERATOR_H
