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

#ifndef VECTOR_FUNCS_H
#define VECTOR_FUNCS_H

#include "gVectorBase.h"

GLINEAR_NAMESPACE_BEGIN

template<typename T>
T append(const T& v1,
         const T& v2);


template<typename T>
void lengthen(T& vec,
              int newlen);

template<typename T>
int findSeq(const T& vec,
            const T& wanted,
            int startAt = 0);
/**
 * Utility used to determine whether a Vector contains nan's. The test is based
 * on the fact that a nan value is not equal to itself (nan != nan). There is no
 * standard method to test for nan's in C++ partly because the representation of
 * a nan value is closely linked to hardware and therefore not easily portable.
 * The test used here is one that is accepted to be correct for most platforms.
 * @param vec gLinear Vector object
 * Return true if Vector contains at least one nan value
 */
template<typename T>
bool hasNan(const gVectorBase<T>& vec);

GLINEAR_NAMESPACE_END

#include "vector_funcs.tcc"

#endif // VECTOR_FUNCS_H
