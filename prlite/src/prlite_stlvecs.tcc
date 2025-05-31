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
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : pre-2000
 * Copyright  : University of Stellenbosch, all rights retained
 */

namespace prlite{

template<typename T>
void setVals(std::vector<T>& vec, const T& first, const T& incr) {
  typename std::vector<T>::size_type sz = vec.size();
  if (sz) {
    vec[0] = first;
    for (decltype(sz) i = 1; i < sz; i++) {
      vec[i] = vec[i-1]+incr;
    } // for i
  } // if
} // set

template<typename T>
std::vector<T>  subrange(const std::vector<T>& vec, unsigned from, unsigned upto) {
  std::vector<T> subvec(vec.begin()+from, vec.begin()+upto+1);
  return subvec;
} // subrange

// find elem
template<typename T>
int findElem(const std::vector<T>& vec,
             const T& wanted,
             int   startAt,
             int direction) {
  // forward
  if (direction >= 0) {
    for (unsigned indx = startAt; indx < vec.size(); indx++) {
      if (vec[indx] == wanted) {
        return indx;
      } // if
    } // for
  } // if
  else { // reverse
    // OK, this is messy. But think it through while you start hacking away.
    // TIP: an unsigned can never be < 0.
    for (int indx = startAt+1; indx-- > 0; ) {
      if (vec[indx] == wanted) {
        return indx;
      } // if
    } // for
  } // else
  // not found
  return -1;
} // findElem

} // namespace prlite

template<typename T>
std::ostream& operator<<(std::ostream& file,
                         const std::vector<T>& vec) {
  typename std::vector<T>::size_type sz = vec.size();
  file << sz << "  ";
  for (unsigned i = 0; i < sz; i++) {
    file << vec[i] << " ";
  } // for i
  return file;
} // operator<<

template<typename T>
std::istream& operator>>(std::istream& file,
                         std::vector<T>& vec) {
  unsigned sz;
  file >> sz;
  vec.resize(sz);
  for (unsigned i = 0; i < sz; i++) {
    file >>  vec[i];
  } // for i
  return file;
} // operator>>

template<typename T>
std::vector<bool> operator<(std::vector<T>& vec, const T& thresh) {
  typename std::vector<T>::size_type sz = vec.size();
  std::vector<bool> answ(sz);
  std::transform(vec.begin(), vec.end(),
                 answ.begin(), [thresh](const T& val) { return val < thresh; } );
  return answ;
} // operator<
