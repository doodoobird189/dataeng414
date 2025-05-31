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
  $Id$
*/

// standard headers
#include <iostream>  // cout, endl

GLINEAR_NAMESPACE_BEGIN

template<typename T_Value>
const char gColVector<T_Value, DENSE>::className[] = "gColVector<T_Value, DENSE>";

template<typename T_Value>
const char gColVector<T_Value, SLICE>::className[] = "gColVector<T_Value, SLICE>";

// No init ctor.
template<typename TF_Value>
inline gColVector<TF_Value, DENSE>::gColVector(gLinearNoInit)
    : T_Base( gLinearNoInit() ) {
} // constructor

// Default ctor and ctor specifying size.
template<typename TF_Value>
inline gColVector<TF_Value, DENSE>::gColVector(T_Index elem)
    : T_Base(elem) {
} // constructor

// Copy ctor.
template<typename TF_Value>
inline gColVector<TF_Value, DENSE>::gColVector(const gColVector<TF_Value,
                                               DENSE>& vec)
    : T_Base(vec) {
} // constructor

// Initializing ctor.
template<typename TF_Value>
inline gColVector<TF_Value, DENSE>::gColVector(T_Index new_nElem,
                                               const T_StorageElem* new_stor,
                                               T_ConstPointer new_startAddr)
    : T_Base(new_nElem, new_stor, new_startAddr) {
} // constructor

template<typename TF_Value>
inline gColVector<TF_Value, DENSE>::gColVector(T_ConstPointer addr,
                 T_Index elem) :
  T_Base(addr, elem) {
}

// Copy assignment. Needed because the operator= doesn't inherit and
// we don't want the compiler to create an implicit one.
template<typename TF_Value>
inline typename gColVector<TF_Value, DENSE>::T_Vector&
gColVector<TF_Value, DENSE>::operator=(const T_Vector& vec) {
  return this->assignTo(vec);
} // operator=

// gColVector assignment. Needed because the operator= doesn't inherit.
template<typename TF_Value>
template<typename T_SrcValue, typename T_SrcStorage>
inline typename gColVector<TF_Value, DENSE>::T_Vector&
gColVector<TF_Value, DENSE>::operator=(const gColVector<T_SrcValue,
               T_SrcStorage>& src) {
  return this->assignTo(src);
} // operator=

template<typename TF_Value>
inline void gColVector<TF_Value, DENSE>::specialInfo(int i) const {
  std::cout << "gColVector<T, DENSE>" << std::endl;
  T_Base::specialInfo(i);
} // specialInfo


/**
 * No init ctor.
 */
template<typename TF_Value>
inline gColVector<TF_Value, SLICE>::gColVector(gLinearNoInit)
    : T_Base( gLinearNoInit() ) {
} // constructor

/**
 * Default ctor and ctor specifying size.
 */
template<typename TF_Value>
inline gColVector<TF_Value, SLICE>::gColVector(T_Index elem)
    : T_Base(elem) {
} // constructor

/**
 * Copy ctor.
 */
template<typename TF_Value>
inline gColVector<TF_Value, SLICE>::gColVector(const gColVector<TF_Value,
                                               SLICE>& vec)
    : T_Base(vec) {
} // constructor

/**
 * Initializing ctor.
 */
template<typename TF_Value>
inline gColVector<TF_Value, SLICE>::gColVector(T_Index new_nElem,
                                               const T_StorageElem* new_stor,
                                               T_ConstPointer new_startAddr,
                                               T_Index new_stride)
    : T_Base(new_nElem, new_stor, new_startAddr, new_stride) {
} // constructor

/**
 * Copy assignment.
 *
 * Needed because the operator= doesn't inherit and we don't want the
 * compiler to create an implicit one.
 */
template<typename TF_Value>
inline typename gColVector<TF_Value, SLICE>::T_Vector&
gColVector<TF_Value, SLICE>::operator=(const T_Vector& vec) {
  return this->assignTo(vec);
} // operator=

/** gColVector assignment. Needed because the operator= doesn't inherit.
 */
template<typename TF_Value>
template<typename T_SrcValue, typename T_SrcStorage>
inline typename gColVector<TF_Value, SLICE>::T_Vector&
gColVector<TF_Value, SLICE>::operator=(const gColVector<T_SrcValue,
               T_SrcStorage>& src) {
  return this->assignTo(src);
} // operator=

template<typename TF_Value>
inline void gColVector<TF_Value, SLICE>::specialInfo(int i) const {
  std::cout << "gColVector<T, SLICE>" << std::endl;
  T_Base::specialInfo(i);
} // specialInfo

GLINEAR_NAMESPACE_END
