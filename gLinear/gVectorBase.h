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

#ifndef GVECTOR_BASE_H
#define GVECTOR_BASE_H

// glinear headers
#include "gLinearMacros.h"
#include "gIndexRange.h"
#include "gVectorTraits.h"
#include "gLinearMacros.h"

// standard headers
#include <iosfwd>  // ostream, istream
#include <string>  // string


GLINEAR_NAMESPACE_BEGIN

/* ________________________________________________________________ */

/**
 * Vector base class defining the interface and some generic algorithms.
 */
template<typename TF_Vector>
class gVectorBase {
public:
  /**
   * Type of vector
   */
  typedef TF_Vector T_Vector;
  /**
   * Type of vector elements.
   */
  typedef typename gVectorTraits<T_Vector>::T_Value T_Value;
  typedef typename gVectorTraits<T_Vector>::T_Value value_type;
  /**
   * Type of pointer to const vector elements.
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstPointer T_ConstPointer;
  /**
   * Type of pointer to vector elements.
   */
  typedef typename gVectorTraits<T_Vector>::T_Pointer T_Pointer;
  /**
   * Type of reference to const vector elements.
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstReference T_ConstReference;
  /**
   * Type of reference to vector elements.
   */
  typedef typename gVectorTraits<T_Vector>::T_Reference T_Reference;
  /**
   * Type of const slice.
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstSlice T_ConstSlice;
  /**
   * Type of slice.
   */
  typedef typename gVectorTraits<T_Vector>::T_Slice T_Slice;
  /**
   * Type of const subrange.
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstSubrange T_ConstSubrange;
  /**
   * Type of subrange.
   */
  typedef typename gVectorTraits<T_Vector>::T_Subrange T_Subrange;
  /**
   * Index type
   */
  typedef typename gVectorTraits<T_Vector>::T_Index T_Index;
  typedef typename gVectorTraits<T_Vector>::T_Index size_type;
  typedef typename gVectorTraits<T_Vector>::T_Index difference_type;
  /**
   * Const Iterator (immuting iterator) type
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstIterator T_ConstIterator;
  /**
   * Iterator type
   */
  typedef typename gVectorTraits<T_Vector>::T_Iterator T_Iterator;
  /**
   * Dense vector type
   */
  typedef typename gVectorTraits<T_Vector>::T_Dense T_Dense;
  /**
   * Const transpose type
   */
  typedef typename gVectorTraits<T_Vector>::T_ConstTranspose T_ConstTranspose;
  /**
   * Transpose type
   */
  typedef typename gVectorTraits<T_Vector>::T_Transpose T_Transpose;
  /**
   * Debug string
   */
  static const char className[];
  /**
   * Info dump starting at most derived
   */
  void info(int i = 0) const;

  /**
   * Return size of vector. \\
   * @return Size of vector.
   */
  T_Index size(void) const;

  /**
   * Resize vector to desired size. \\
   * The vector is only changed if the requested size {\em new_size} differs
   * from it's current size. \\
   * WARNING: Do not attempt this on a vector containing aliases. It won't work :)
   * @param new_size The requested size.
   * @return
   * \begin{itemize}
   * \item {\bf true -} Resize changed the object.
   * \item {\bf false -} Resize returned without modifying this object.
   * \end{itemize}
   */
  bool resize(T_Index elem);

  /**
   * Return the const transpose of the vector. \\
   */
  inline const T_Transpose transpose(void) const;
  /**
   * Return the transpose of the vector. \\
   */
  inline T_Transpose transpose(void);
protected:
  /**
   * Assignment. \\
   * Performs a deep copy from the argument vector to this vector.
   * First attempts a resize, which would fail if any alises are present
   * on this vector.
   * @param vec Source vector
   * @returns *this
   */
  template<typename T_SrcVector>
  T_Vector& assignTo(const gVectorBase<T_SrcVector>& vec);

public:

  /**
   * Scalar assignment operator. \\
   * Assign the scalar argument to each value of the gVectorBase.
   * @param scalar Scalar argument.
   * @return *this.
   */
  T_Vector& assignToAll(T_ConstReference scalar);

  /**
   * Assign a set of values to this vector. \\
   * Assign the set consisting of {\em nElem} values, each {\em increment}
   * apart, starting at {\em beginval} to this vector.
   * Assignment and addition should be overloaded for class T.
   * @param beginval Start value.
   * @param increment Increment.
   */
  T_Vector& set(T_ConstReference beginval,
                T_ConstReference increment);

  /**
   * Const element access. \\
   */
  T_ConstReference operator[](T_Index elem) const;
  /**
   * Element access. \\
   * Return reference to the elem'th element of the vector. \\
   * @param elem Index of element to access. Valid values 0 to nElem-1.
   * @return Reference to the elem'th element.
   */
  T_Reference operator[](T_Index elem);

  /**
   * Assign values from C/C++ array.
   */
  void assignFromArray(T_ConstPointer addr,
                       T_Index elem = -1,
                       T_Index stride = 0);

  /**
   * Return const iterator referring to the first element. \\
   * This iterator may be used to traverse and inspect
   * the vector.
   */
  T_ConstIterator begin(void) const;
  /**
   * Return iterator referring to the first element. \\
   * This iterator may be used to traverse, inspect and modify
   * the vector.
   */
  T_Iterator begin(void);
  /**
   * Return the past-the-end iterator. \\
   * If an iterator equals this value, it should {\bf NOT} be dereferenced
   * or incremented.
   */
  T_ConstIterator end(void) const;
  /**
   * Return the past-the-end iterator. \\
   * If an iterator equals this value, it should {\bf NOT} be dereferenced
   * or incremented.
   */
  T_Iterator end(void);

  /**
   * Return subrange slice of vector. \\
   * The returned vector is a shallow copy of this vector.
   * @param elems Parameter specifying the start index, end index and
   * stride of the slice.
   * @return T_Slice object addressing a subset of the data that this vector does.
   */
  T_Slice slice(const gIndexRange& elems);
  /**
   * Return const subrange slice of vector. \\
   * The returned vector is a shallow copy of this vector.
   * @param elems Parameter specifying the start index, end index and
   * stride of the slice.
   */
  T_ConstSlice slice(const gIndexRange& elems) const;
  /**
   * Return subrange of vector.
   * @return fromElem Index of first element of subrange.
   * @return toElem Index of last element of subrange.
   */
  T_Subrange subrange(T_Index fromElem,
                      T_Index toElem);
  /**
   * Return const subrange of vector.
   * @return fromElem Index of first element of subrange.
   * @return toElem Index of last element of subrange.
   */
  T_ConstSubrange subrange(T_Index fromElem,
                           T_Index toElem) const;

  /**
   * Return the number of references (aliases) existing on this object.
   */
  T_Index getRefs(void) const;

  /**
   * Ensure that no aliases exist to the data contained in this object.
   */
  void deepen(void) const;

  /**
   * Let this gVectorBase be an alias of another one.
   */
  void alias(T_Vector& vec);

  /**
   * Let this gVectorBase be an alias of another (const) one.
   */
  void alias(const T_Vector& vec) const;

  /* Generic Algorithms based on iterators */

  /**
   * Return value of the minimum (smallest) element. \\
   * The T_Value::operator<(const T_Value& elem) function call is used to find the
   * smallest element.
   * @param start Start index (default: start of vector).
   * @param end End index (default: end of vector).
   * @return The minimum-valued element.
   */
  T_Value min(void) const;
  /*T_Index start = -1, T_Index end = -1*/

  /**
   * Return value of the of the minimum (smallest) element. \\
   * The T_Value::operator<(const T_Value& elem) function call is used to find the
   * smallest element. If more than one element shares the same minimum
   * value, the one with the highest index is returned.
   * @param start Start index (default: start of vector).
   * @param end End index (default: end of vector).
   * @return The minimum-valued element.
   */
  T_Value minVal(T_Index start = 0, T_Index end = -1) const;

  /**
   * Return index of the of the minimum (smallest) element. \\
   * The T_Value::operator<(const T_Value& elem) function call is used to find the
   * smallest element. If more than one element shares the same minimum
   * value, the one with the highest index is returned.
   * @param start Start index (default: start of vector).
   * @param end End index (default: end of vector).
   * @return Highest index corresponding to the minimum-valued element.
   */
  T_Index minPos(T_Index start = 0, T_Index end = -1) const;

  /**
   * Return value of the maximum (largest) element. \\
   * The T_Value::operator>(const T_Value& elem) function call is used to find the
   * largest element.
   * @return The maximum-valued element.
   */
  T_Value max(void) const;

  /**
   * Return value of the of the maximum (largest) element. \\
   * The T_Value::operator>(const T_Value& elem) function call is used to find the
   * largest element. If more than one element shares the same maximum
   * value, the one with the highest index is returned.
   * @param start Start index (default: start of vector).
   * @param end End index (default: end of vector).
   * @return The maximum-valued element.
   */
  T_Value maxVal(T_Index start = 0, T_Index end = -1) const;

  /**
   * Return index of the of the maximum (largest) element. \\
   * The T_Value::operator>(const T_Value& elem) function call is used to find the
   * largest element. If more than one element shares the same maximum
   * value, the one with the highest index is returned.
   * @param start Start index (default: start of vector).
   * @param end End index (default: end of vector).
   * @return Highest index corresponding to the maximum-valued element.
   */
  T_Index maxPos(T_Index start = 0, T_Index end = -1) const;

  /**
   * Return value of the maximum (largest) magnitude element. \\
   * The T_Value::operator>(const T_Value& elem) function call is used to find the
   * largest magnitude.
   * @return The maximum-valued magnitude.
   */
  T_Value absMax(void) const;

  /**
   * Search for element in vector and return index. \\
   * Search through the vector for the first occurence of the argument
   * and return the index. A -1 is returned if the element is not found.
   * @param wanted Element to search for.
   * @param startAt Start index. Default = 0.
   * @param Direction {\bf >=0} Forward (default) \\
   * {\bf <0 } Backward.
   * @return {\bf -1 -} Element not found. \\
   * {\bf n -} Index of element.
   */
  T_Index findElem(T_ConstReference wanted,
                   T_Index startAt = 0,
                   int direction = 1) const;

  /** Vector computed assignment: addition
   */
  template<typename T_SrcVector>
  T_Vector& operator+=(const gVectorBase<T_SrcVector>& vec);
  /** Vector computed assignment: subtraction
   */
  template<typename T_SrcVector>
  T_Vector& operator-=(const gVectorBase<T_SrcVector>& vec);
  /** Scalar computed assignment: addition
   */
  T_Vector& operator+=(T_ConstReference scalar);
  /** Scalar computed assignment: subtraction
   */
  T_Vector& operator-=(T_ConstReference scalar);
  /** Scalar computed assignment: multiplication
   */
  T_Vector& operator*=(T_ConstReference scalar);
  /** Scalar computed assignment: division
   */
  T_Vector& operator/=(T_ConstReference scalar);

  /**
   * Return sum of all elements.
   * Calculate the sum of all the elements using the
   * T_Value::operator+=(const T_Value& t) function call.
   * @return Sum of all elements.
   */
  T_Value sum(void) const;

  /**
   * Return sum of square of elements. \\
   * This is equal to the square of the norm of the vector. \\
   * The square is calculated using T_Value::operator*(const T_Value& t) and
   * the accumulation is performed using T_Value::operator+=(const T_Value& t)
   * @return Sum of square of element.
   */
  T_Value sum2(void) const;

  /**
   * Print contents to cout
   */
  void print(int depth = 0) const;
  /**
   * Binary Output
   * Write the length of the vector, followed by the contets
   * to a stream (binary output).
   */
  std::ostream& write(std::ostream& os) const;
  /**
   * Binary Input
   * Read the length of the vector, followed by the contets
   * from a stream (binary input).
   */
  std::istream& read(std::istream& is);
private: // support
  void chkidxIndex(std::string funcName, std::string infostr, T_Index elem) const;
  void chkidxZerosize(std::string funcName, std::string infostr) const;
  void chkidxSize(std::string funcName, std::string infostr, T_Index size1, T_Index size2) const;
  void chkidxSubrange(std::string funcName, std::string infostr, T_Index from, T_Index to) const;
  void chkidxSlice(std::string funcName, std::string infostr, const gIndexRange& elems) const;
  void chkidxAlias(std::string funcName, std::string infostr) const;
}; // gVectorBase

GLINEAR_NAMESPACE_END

#include "gVectorBase.tcc"

#endif // GVECTOR_BASE_H
