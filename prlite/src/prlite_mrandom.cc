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

char MRANDOM_CC_VERSION[] = "$Id$";

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14/09/95
         PURPOSE: Random number generator.

*******************************************************************************/

// patrec headers
#include "prlite_mrandom.hpp"
#include "prlite_stdfun.hpp"  // M_PI
#include "prlite_logging.hpp"

// standard headers
#include <iostream>
#include <cmath>  // cos, log, sqrt

using namespace std;

/******************************************************************************/

namespace prlite{

MRandom offset(13);

MRandom::MRandom() {
  randomize();
} // default constructor

MRandom::MRandom(long x) {
  if (x < 1 || x > 2147483647L) {
    PRLITE_CAUTION("random::random : invalid seed" << endl);
  } // if
  seed = x;
  Flat();
} // constructor

void MRandom::randomize() {
  time(&seed);
  seed += static_cast<long>( offset.Int(32000) );
  if (seed <= 0) {
    seed += 2147483647L;
  } // if
  Flat();
} // randomize

float MRandom::Flat() {
  static const long int a = 16807L;          // multiplier
  static const long int m = 2147483647L;     // modulo
  static const long int q = 127773L;         // m div a
  static const long int r = 2836L;           // m mod a
  long test;
  ldiv_t z;
  z = ldiv(seed, q);
  test = a * z.rem - r * z.quot;
  if (test > 0) {
    seed = test;
  } // if
  else {
    seed = test + m;
  } // else
  return static_cast<float>(seed) / static_cast<float>(m);
} // Flat

float MRandom::Flat(float from,
                    float upto) {
  return Flat() * (upto - from) + from;
} // Flat

int MRandom::Int(int max) {
  if (max < 0) {
    PRLITE_CAUTION("random::Int : maxVal < 0" );
  } // if
  int tmp = static_cast<int>(Flat() * static_cast<float>(max) + 0.5); // rounded
  if (tmp == max) {
    tmp = 0;
  } // if
  return (tmp);
} // Int

int MRandom::Int(int from,
                 int upto) {
  if (upto < from) {
    PRLITE_CAUTION("random::Int : [" << from << ':' << upto << "]");
  } // if
  int tmp = static_cast<int>( Flat((float) from, (float) upto+1) ); // rounded
  if (tmp > upto) {
    tmp = from;
  } // if
  return (tmp);
} // Int

float MRandom::Gauss() {
  return ( sqrt( -2.0f * log( Flat() ) ) * cos( 2.0f * static_cast<float>(M_PI) * Flat() ) );
} // Gauss

float MRandom::Gauss(float mean,
                     float stdev) {
  if (stdev < 0) {
    PRLITE_CAUTION("random::Gauss : stdev < 0.0");
  } // if
  return ( mean + stdev * Gauss() );
} // Gauss

} // namespace prlite
