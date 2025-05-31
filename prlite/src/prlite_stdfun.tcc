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

// patrec headers
#include "prlite_error.hpp"  // FAIL

// standard headers
#include <limits>  // numeric_limits
#include <cstdlib>  // abs
#include <algorithm>  // min

namespace prlite{

template<typename T>
int round2int(T x) {
  if ( x < -std::numeric_limits<int>::max() || x > std::numeric_limits<int>::max() ) {
    PRLITE_FAIL("round2int: " << x << " out of range");
  } // if
  return x >= 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
} // round2int

template<typename T>
long round2long(T x) {
  if ( x < -std::numeric_limits<long int>::max() || x > std::numeric_limits<long int>::max() ) {
    PRLITE_FAIL("round2long: " << x << " out of range");
  } // if
  return x >= 0 ? static_cast<long>(x + 0.5) : static_cast<long>(x - 0.5);
} // round2long

template<typename T>
unsigned round2unsigned(T x) {
  if ( x < 0 || x > std::numeric_limits<unsigned int>::max() ) {
    PRLITE_FAIL("round2uint: " << x << " out of range");
  } // if
  return static_cast<unsigned>(x + 0.5);
} // round2unsigned

template<typename T>
unsigned long round2ulong(T x) {
  if ( x < 0 || x > std::numeric_limits<unsigned long int>::max() ) {
    PRLITE_FAIL("round2ulong: " << x << " out of range");
  } // if
  return static_cast<unsigned long>(x + 0.5);
} // round2ulong

template<typename T>
bool realEqual(T x,
               T y) {
  T a( std::abs(x) );
  T b( std::abs(y) );
  T m( std::min(a, b) );
  return std::abs(x - y) <= (std::numeric_limits<T>::epsilon() * m);
} // realEqual

// ----------------------------------------------------------------------------
// Combination : returns a combination b
template<typename T>
T combination(int a, int b) {
  T tp = static_cast<T>(a);
  T bp = static_cast<T>(1.0);

  if (b == 0) {
    return 1;
  } // if

  if (a == 0) {
    return 0;
  } // if

  for (int i = 1; i < b; i++) {
    tp *= static_cast<T>(a - i);
    bp *= static_cast<T>(i + 1);
  } // for

  return (tp / bp);
} // combination

} // namespace prlite
