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

GLINEAR_NAMESPACE_BEGIN

// The following two functions (transformInto) are adapted from stl_algo.h,
// included in the STL release shipped with egcs-1.1.2-release (egcs-2.91.66 19990314)
// Copyright retained and displayed as required by license:
// Copyright (c) 1994
// Hewlett-Packard Company
//
// Copyright (c) 1996
// Silicon Graphics Computer Systems, Inc.

// transformInto: unary transfrom
// calling binary (void-return, unary arg) functor.
template<typename T_InputIter, typename T_OutputIter, typename T_UnaryIntoFunctor>
T_OutputIter transformInto(T_InputIter first,
                           T_InputIter last,
                           T_OutputIter result,
                           T_UnaryIntoFunctor func) {
  for ( ; first != last; ++first, ++result) {
    func(*first, *result);
  } // for
  return result;
} // transformInto

// transformInto: binary transform
// calling ternary (void-return, binary arg) functor.
template<typename T_InputIter1, typename T_InputIter2, typename T_OutputIter, typename T_BinaryIntoFunctor>
T_OutputIter transformInto(T_InputIter1 first1,
                           T_InputIter1 last1,
                           T_InputIter2 first2,
                           T_OutputIter result,
                           T_BinaryIntoFunctor func) {
  for ( ; first1 != last1; ++first1, ++first2, ++result) {
    func(*first1, *first2, *result);
  } // for
  return result;
} // transformInto

GLINEAR_NAMESPACE_END
