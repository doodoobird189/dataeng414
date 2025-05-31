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
  DEPRECATED! rather use the new C++17 functions
*/

#ifndef PRLITE_MRANDOM_HPP
#define PRLITE_MRANDOM_HPP

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14/09/95
         PURPOSE: Random number generator.

*******************************************************************************/

#include <ctime>  // time

namespace prlite{
/**
  The class MRandom supplies the routines for supplying basic random number
  generation.
*/

class MRandom {

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
public:

  /// randomised seed
  MRandom();
  /// using seed x
  MRandom(long x);

  //========================= Operators and Conversions =========================

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

  //=========================== Inplace Configuration ===========================

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================

  //======================= Non-virtual Member Functions ========================
public:

  /// system clock in random seed
  void randomize();

  /// flat distrib (0:1]
  float Flat();

  /// flat distrib between limits
  float Flat(float from,
             float upto);

  /// discr distrib [0:max-1]
  int Int(int max);

  /// discr distrib [from:upto]
  int Int(int from,
          int upto);

  /// Gaussian distrib 0 mean, stdev 1
  float Gauss();

  /// Gaussian distrib
  float Gauss(float mean,
              float stdev);

  //================================== Friends ==================================

  //=============================== Data Members ================================
private:

   time_t seed;                            // any between 1 and 2147483647

  //============================ DEPRECATED Members =============================


}; // MRandom

} // namespace prlite

#endif // PRLITE_MRANDOM_HPP
