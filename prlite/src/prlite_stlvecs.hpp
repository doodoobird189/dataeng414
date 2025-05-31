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

/*******************************************************************************
      AUTHOR: J.A du Preez (Copyright University of Stellenbosch, all rights reserved.)

Purpose:
This code supplies a number of functions/operaters on an stl::vector
to broadly duplicate some capability available in the gLinear vectors.
*******************************************************************************/

#ifndef PRLITE_STLVECS_HPP
#define PRLITE_STLVECS_HPP

// standard headers
#include <vector>
#include <algorithm>

 /******************************************************************************/

namespace prlite{

 /**
  * Initialises the vector starting with 'first' in vec[0] and
  * thereafter incrementing the value with 'incr' in each subsequent
  * element. (gLinear: vec.set)
  */
template<typename T>
void setVals(std::vector<T>& vec, const T& first, const T& incr);

 /**
  * Returns the index of the first occurence of 'wanted', or -1 if not
  * found. (gLinear: vec.findElem)
  */
template<typename T>
int findElem(const std::vector<T>& vec,
             const T& wanted,
             int   startAt=0,
             int direction=1);

/**
 * Returns a vector consisting of the elements between (and including)
 * 'from' to 'upto'.  NOTE, upto is included in the range (i.e. it is not
 * past-end). CAREFUL, different from gLinear's subrange this is not a
 * reference/slice on the underlying array. (gLinear: vec.subrange)
 */
template<typename T>
std::vector<T> subrange(const std::vector<T>& vec, unsigned from, unsigned upto);


} // namespace prlite

/**
 * writes vector to ostream. (gLinear: <<)
 */
template<typename T>
std::ostream& operator<<(std::ostream& file,
                         const std::vector<T>& vec);
/**
 * reads vector from istream. (gLinear: >>)
 */
template<typename T>
std::istream& operator>>(std::istream& file,
                         std::vector<T>& vec);

/**
 * returns a boolean vector of same size as 'vec' with 0/false where
 * condition is false and 1/true where condition is true.
 */
template<typename T>
std::vector<bool> operator<(std::vector<T>& vec, const T& thresh);

/******************************************************************************/

#include "prlite_stlvecs.tcc"


#endif // PRLITE_STLVECS_HPP
