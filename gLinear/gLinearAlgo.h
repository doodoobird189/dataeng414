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

#ifndef GLINEAR_ALGO_H
#define GLINEAR_ALGO_H

GLINEAR_NAMESPACE_BEGIN

// transform into: unary func
template<typename T_InputIter, typename T_OutputIter, typename T_UnaryIntoFunctor>
T_OutputIter transformInto(T_InputIter first1,
                           T_InputIter last1,
                           T_OutputIter result,
                           T_UnaryIntoFunctor func);

// transform into: binary func
template<typename T_InputIter1, typename T_InputIter2, typename T_OutputIter, typename T_BinaryIntoFunctor>
T_OutputIter transformInto(T_InputIter1 first1,
                           T_InputIter1 last1,
                           T_InputIter2 first2,
                           T_OutputIter result,
                           T_BinaryIntoFunctor func);

GLINEAR_NAMESPACE_END

#include "gLinearAlgo.tcc"

#endif // GLINEAR_ALGO_H
