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
const char gRowVector<T_Value, DENSE>::className[] = "gRowVector<T_Value, DENSE>";



template<typename T_Value>
const char gRowVector<T_Value, SLICE>::className[] = "gRowVector<T_Value, SLICE>";



template<typename TF_Value>
inline gRowVector<TF_Value, DENSE>::gRowVector(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, DENSE>::gRowVector(T_Index elem) :
  T_Base(elem) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, DENSE>::gRowVector(const gRowVector<TF_Value, DENSE>& vec) :
  T_Base(vec) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, DENSE>::gRowVector(T_Index new_nElem,
                                               const T_StorageElem* new_stor,
                                               T_ConstPointer new_startAddr) :
  T_Base(new_nElem, new_stor, new_startAddr) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, DENSE>::gRowVector(T_ConstPointer addr,
                 T_Index elem) :
  T_Base(addr, elem) {
}

template<typename TF_Value>
inline typename gRowVector<TF_Value, DENSE>::T_Vector&
gRowVector<TF_Value, DENSE>::operator=(const T_Vector& vec) {
  return this->assignTo(vec);
}

template<typename TF_Value>
template<typename T_SrcValue, typename T_SrcStorage>
inline typename gRowVector<TF_Value, DENSE>::T_Vector&
gRowVector<TF_Value, DENSE>::operator=(const gRowVector<T_SrcValue, T_SrcStorage>& src) {
  return this->assignTo(src);
}

template<typename TF_Value>
inline void gRowVector<TF_Value, DENSE>::specialInfo(int i) const {
  std::cout << "gRowVector<T, DENSE>" << std::endl;
  T_Base::specialInfo(i);
}

template<typename TF_Value>
inline gRowVector<TF_Value, SLICE>::gRowVector(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, SLICE>::gRowVector(T_Index elem) :
  T_Base(elem) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, SLICE>::gRowVector(const gRowVector<TF_Value, SLICE>& vec) :
  T_Base(vec) {
}

template<typename TF_Value>
inline gRowVector<TF_Value, SLICE>::gRowVector(T_Index new_nElem,
                                               const T_StorageElem* new_stor,
                                               T_ConstPointer new_startAddr,
                                               T_Index new_stride) :
  T_Base(new_nElem, new_stor, new_startAddr, new_stride) {
}

template<typename TF_Value>
inline typename gRowVector<TF_Value, SLICE>::T_Vector&
gRowVector<TF_Value, SLICE>::operator=(const T_Vector& vec) {
  return this->assignTo(vec);
}

template<typename TF_Value>
template<typename T_SrcValue, typename T_SrcStorage>
inline typename gRowVector<TF_Value, SLICE>::T_Vector&
gRowVector<TF_Value, SLICE>::operator=(const gRowVector<T_SrcValue, T_SrcStorage>& src) {
  return this->assignTo(src);
}

template<typename TF_Value>
inline void gRowVector<TF_Value, SLICE>::specialInfo(int i) const {
  std::cout << "gRowVector<T, SLICE>" << std::endl;
  T_Base::specialInfo(i);
}

GLINEAR_NAMESPACE_END
