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

/*
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : pre-2000
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_genvec.hpp"  // Vector

namespace prlite{

template<typename T>
T append( const T& v1,
    const T& v2 ) {
  T tmp( v1.size() + v2.size() );
  int idx(0);
  decltype(v1.size()) i;
  for (i = 0; i < v1.size(); i++) {
    tmp[idx++] = v1[i];
  } // for i
  for (i = 0; i < v2.size(); i++) {
    tmp[idx++] = v2[i];
  } // for i
  return tmp;
} // append


template<typename T>
void lengthen( T& vec,
         int newlen ) {
  T tmp;
  tmp = vec;
  vec.resize(newlen);
  for (int i = 0; i < tmp.size(); i++) {
    if ( i == vec.size() ) {
      break;
    } // if
    vec[i] = tmp[i];
  } // for i
} // lengthen


template<typename T>
int findSeq( const T& vec,
       const T& wanted,
       int startAt ) {
  int pos;
  do {
    pos = vec.findElem(wanted[0], startAt);
    startAt = pos + 1;
    if ( pos < 0 || pos + wanted.size() > vec.size() ) {
      return -1;
    } // if
    int i;
    for (i = 1; i < wanted.size(); i++) {
      if ( vec[pos + i] != wanted[i] ) {
        break;
      } // if
    } // for i
    if ( i == wanted.size() ) {
      return pos;
    } // if
  } while (pos >= 0);
  return -1;
} // findSeq


template<typename VEC>
VEC extract( const VEC& vec,
             const std::vector<int>& idxs ) {

  VEC ret( idxs.size() );
  for (unsigned i = 0; i < static_cast<unsigned>( idxs.size() ); i++) {
    ret[i] = vec[ idxs[i] ];
  } // for i

  return ret;

} // extract

template<typename VEC>
VEC extract( const VEC& vec,
             const std::vector<unsigned>& idxs ) {

  VEC ret( idxs.size() );
  for (unsigned i = 0; i < static_cast<unsigned>( idxs.size() ); i++) {
    ret[i] = vec[ idxs[i] ];
  } // for i

  return ret;

} // extract

template<typename VEC>
std::vector<int> find( const VEC& vec) {

  unsigned sz = vec.size();
  std::vector<int> ret;
  ret.reserve(sz);
  for (unsigned i = 0; i < sz; i++) {
      if (vec[i]) {
        ret.push_back(i);
      } // if
  } // for i

  return ret;

} // find

  template<typename T>
  gLinear::gRowVector<T> asVector(const std::vector<T>& stlVec) {
    int sz = stlVec.size();
    gLinear::gRowVector<T> prVec(sz);

    const T* stlPtr = &stlVec[0];
    T* prPtr = &prVec[0];
    while (sz--) {
      *prPtr++ = *stlPtr++;
    } // while

    return prVec;
  } // asVector


template<typename T>
gLinear::gRowVector<T> vector2Vector(const std::vector<T>& stlVec) {
  int sz = stlVec.size();
  gLinear::gRowVector<T> prVec(sz);

  const T* stlPtr = &stlVec[0];
  T* prPtr = &prVec[0];
  while (sz--) {
    *prPtr++ = *stlPtr++;
  } // while

  return prVec;
} // vector2Vector

template<typename T>
gLinear::gColVector<T> vector2ColVector(const std::vector<T>& stlVec) {
  int sz = stlVec.size();
  gLinear::gColVector<T> prVec(sz);

  const T* stlPtr = &stlVec[0];
  T* prPtr = &prVec[0];
  while (sz--) {
    *prPtr++ = *stlPtr++;
  } // while

  return prVec;
} // vector2ColVector

template<typename T>
std::vector<T> Vector2vector(const gLinear::gRowVector<T>& glVec) {
  unsigned int sz = glVec.size();
  std::vector<T> stlVec;

  for (unsigned int i(0); i < sz; i++) {
    stlVec.push_back( glVec[i] );
  } // for

  return stlVec;
} // Vector2vector

template<typename T>
std::vector<T> ColVector2vector(const gLinear::gColVector<T>& glVec) {
  unsigned int sz = glVec.size();
  std::vector<T> stlVec;

  for (unsigned int i(0); i < sz; i++) {
    stlVec.push_back( glVec[i] );
  } // for

  return stlVec;
} // ColVector2vector

} // namespace prlite
