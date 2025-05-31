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

template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(gLinearNoInit) :
  T_Base( gLinearNoInit() ) {
}

template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(void) :
  T_Base() {
}


template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(T_Index _rows,
                                              T_Index _cols) :
  T_Base(_rows, _cols) {
}

template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(const T_Matrix& mat) :
  T_Base(mat) {
}

template<typename TF_Value, typename TF_Storage>
template<typename T_VecStor>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(const gRowVector<T_Value, T_VecStor>& vec) :
  T_Base(vec) {

}

template<typename TF_Value, typename TF_Storage>
template<typename T_VecStor>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(const gColVector<T_Value, T_VecStor>& vec) :
  T_Base(vec) {
}

template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>::gMatrix(T_ConstPointer addr,
                                              T_Index rows, T_Index cols) :
  T_Base(addr, rows, cols) {
}

template<typename TF_Value, typename TF_Storage>
inline gMatrix<TF_Value, TF_Storage>&
gMatrix<TF_Value, TF_Storage>::operator=(const gMatrix& mat) {
  return this->assignTo(mat);
}

template<typename TF_Value, typename TF_Storage>
template<typename T_SrcValue, typename T_SrcStorage>
inline gMatrix<TF_Value, TF_Storage>&
gMatrix<TF_Value, TF_Storage>::operator=(const gMatrix<T_SrcValue, T_SrcStorage>& mat) {
  return this->assignTo(mat);
}

template<typename TF_Value, typename TF_Storage>
inline void gMatrix<TF_Value, TF_Storage>::specialInfo(int i) const {
  std::cout << "gMatrix<T_Value, " << T_Storage::className << ">" << std::endl;
}

GLINEAR_NAMESPACE_END
