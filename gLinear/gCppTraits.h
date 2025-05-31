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

#ifndef GCPP_TRAITS_H
#define GCPP_TRAITS_H

// standard headers
#include <complex>  // complex

/**
 * Traits relating standard C++ types
 */
template<typename T>
struct gCppTraits;

template<>
struct gCppTraits<float> {
  enum {
    isReal = true,
    isComplex = false
  };
  typedef float T_Real;
  typedef std::complex<float> T_Complex;
}; // gCppTraits

template<>
struct gCppTraits< std::complex<float> > {
  enum {
    isReal = false,
    isComplex = true
  };
  typedef float T_Real;
  typedef std::complex<float> T_Complex;
}; // gCppTraits

template<>
struct gCppTraits<double> {
  enum {
    isReal = true,
    isComplex = false
  };
  typedef double T_Real;
  typedef std::complex<double> T_Complex;
}; // gCppTraits

template<>
struct gCppTraits< std::complex<double> > {
  enum {
    isReal = true,
    isComplex = false
  };
  typedef double T_Real;
  typedef std::complex<double> T_Complex;
}; // gCppTraits



/**
 * C++ Precision Traits classes
 */

template<typename T1, typename T2>
class gNumPrec2 {
public:
  typedef T1 T_Highest;  // first arg. default highest
  typedef T2 T_Lowest;  // second arg. default lowest
}; // gNumPrec2

// T, T
template<typename T>
class gNumPrec2<T, T> {
public:
  typedef T T_Highest;
  typedef T T_Lowest;
}; // gNumPrec2

// float, double
template<>
class gNumPrec2<float, double> {
public:
  typedef double T_Highest;
  typedef float T_Lowest;
}; // gNumPrec2

// double, float
template<>
class gNumPrec2<double, float> {
public:
  typedef double T_Highest;
  typedef float T_Lowest;
}; // gNumPrec2

#endif // GCPP_TRAITS_H
