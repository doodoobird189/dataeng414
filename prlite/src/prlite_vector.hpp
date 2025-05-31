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

/*******************************************************************************
      AUTHOR: J.A du Preez (Copyright University of Stellenbosch, all rights reserved.)
*******************************************************************************/

#ifndef PRLITE_VECTOR_HPP
#define PRLITE_VECTOR_HPP

// patrec headers
#include "prlite_genvec.hpp"  // Vector

// standard headers
#include <vector>
/******************************************************************************/

namespace prlite{

///
template<typename T>
T append(const T& v1,
         const T& v2);

///
template<typename T>
void lengthen(T& vec,
              int newlen);

///
template<typename T>
int findSeq(const T& vec,
            const T& wanted,
            int startAt = 0);

/// A temporary hack until gLinear supports a gIndexedVector
template<typename VEC>
VEC extract(const VEC& vec,
            const std::vector<int>& idxs);
template<typename VEC>
VEC extract(const VEC& vec,
            const std::vector<unsigned>& idxs);

template<typename VEC>
std::vector<int> find( const VEC& vec);

/// converts a std::vector to a gLinear Vector
template<typename T>
gLinear::gRowVector<T> asVector(const std::vector<T>& stlVec);

template<typename T>
gLinear::gRowVector<T> vector2Vector(const std::vector<T>& stlVec);

template<typename T>
gLinear::gColVector<T> vector2ColVector(const std::vector<T>& stlVec);

template<typename T>
std::vector<T> Vector2vector(const gLinear::gRowVector<T>& glVec);

template<typename T>
std::vector<T> ColVector2vector(const gLinear::gColVector<T>& glVec);

/******************************************************************************/

} // namespace prlite

#include "prlite_vector.tcc"

#endif // PRLITE_VECTOR_HPP
