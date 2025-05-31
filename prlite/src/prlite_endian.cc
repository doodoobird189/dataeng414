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

char ENDIAN_CC_VERSION[] = "$Id$";

 /*****************************************************************************
 *               |                                                            *
 * Author        | Dr J.A. du Preez (typically)                               *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 2002/03/07                                                 *
 * Description   | Swap bytes to convert endianness                           *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

// patrec headers
#include "prlite_endian.hpp"

// standard headers
#include <algorithm>

using namespace std;

namespace prlite{

/********************* Endian swapping routines ******************************/

const bool BIG_ENDIAN_CPU(!isLittleEndianCPU() );

bool isLittleEndianCPU() {
  short test = 1;
  char* byte0( reinterpret_cast<char *>(&test) );
  return *byte0 != 0;
} // isLittleEndian

// and some specific ones for efficiency

float endianSwap(float& val) {
  float a, b;
  a = 1;
  b = 2;
  swap(a, b);
  char* byteArr( reinterpret_cast<char *>(&val) );
  swap( byteArr[0], byteArr[3] );
  swap( byteArr[1], byteArr[2] );
  return val;
} // endianSwap

int endianSwap(int& val) {
  char* byteArr( reinterpret_cast<char *>(&val) );
  swap( byteArr[0], byteArr[3] );
  swap( byteArr[1], byteArr[2] );
  return val;
} // endianSwap

unsigned endianSwap(unsigned& val) {
  char* byteArr( reinterpret_cast<char *>(&val) );
  swap( byteArr[0], byteArr[3] );
  swap( byteArr[1], byteArr[2] );
  return val;
} // endianSwap

short endianSwap(short& val) {
  char* byteArr( reinterpret_cast<char *>(&val) );
  swap( byteArr[0], byteArr[1] );
  return val;
} // endianSwap

// bool isLittleEndian();

} // namespace prlite
