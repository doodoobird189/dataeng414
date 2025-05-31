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

GLINEAR_NAMESPACE_BEGIN

// from J.A. du Preez's original patrec code.
template<typename T>
T append(const T& v1,
         const T& v2) {
  T tmp( v1.size() + v2.size() );
  int idx = 0;
  int i;
  for (i = 0; i < v1.size(); i++) {
    tmp[idx++] = v1[i];
  } // for i
  for (i = 0; i < v2.size(); i++) {
    tmp[idx++] = v2[i];
  } // for i
  return tmp;
} // append


template<typename T>
void lengthen(T& vec,
              int newlen) {
  T tmp(vec); // shallow
  vec.deepen(); // branch off
  vec.resize(newlen); // resize and loose data
  for (int i = 0; i < tmp.size(); i++) {
    if ( i == vec.size() ) {
      break;
    } // if
    vec[i] = tmp[i];
  } // for i
} // lengthen

// from J.A. du Preez's original patrec code.
template<typename T>
int findSeq(const T& vec,
            const T& wanted,
            int startAt) {
  int pos;
  do {
    pos = vec.findElem(wanted[0], startAt);
    startAt = pos + 1;
    if ( pos < 0 || pos + wanted.size() > vec.size() ) {
      return -1;
    } // if
    int i;
    for (i = 1; i < wanted.size(); i++) {
      if ( vec[pos+i] != wanted[i] ) {
        break;
      } // if
    } // for i
    if ( i == wanted.size() ) {
      return pos;
    } // if
  } while (pos >= 0);
  return -1;
} // findSeq

template<typename T>
bool hasNan(const gVectorBase<T>& vec) {
  typedef typename gVectorBase<T>::T_Index T_Index;
  for (T_Index i=0; i< vec.size(); i++) {
     if (vec[i] != vec[i]) //test for nan value
     return true;
  } // for
  return false;
} // isValid

GLINEAR_NAMESPACE_END
