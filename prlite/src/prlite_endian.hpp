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

#ifndef PRLITE_ENDIAN_HPP
#define PRLITE_ENDIAN_HPP

namespace prlite{

/********************* Endian swapping routines ******************************/


/*
 * Converts an input parameter in the endian-ness of the current cpu to little endian.
 *
 * @param val a parameter with the endian-ness of the current cpu archictecture.
 * @return the little endian version of val.
 */
template<typename T>
T toLittleEndian(const T& val);

/*
 * Converts an input parameter in the endian-ness of the current cpu to big endian.
 *
 * @param val a parameter with the endian-ness of the current cpu archictecture.
 * @return the big endian version of val.
 */
template<typename T>
T toBigEndian(const T& val);

/*
 * Converts an input parameter in little endian to the endian-ness of the current cpu.
 *
 * @param val a parameter in little endian format.
 * @return val in th e the endian-ness of the current cpu archictecture.
 */
#define PRLITE_fromLittleEndian prlite::toLittleEndian

/*
 * Converts an input parameter in big endian to the endian-ness of the current cpu.
 *
 * @param val a parameter in big endian format.
 * @return val in th e the endian-ness of the current cpu archictecture.
 */
#define PRLITE_fromBigEndian toBigEndian

/**
 * constant is true if machine is big endian, false otherwise
 * would be faster to use an ifdef on the architecture!
 */
extern const bool BIG_ENDIAN_CPU;

/**
 * Tests the endian-ness of the cpu.
 * @return true if machine is little endian, false otherwise.
 */
extern bool isLittleEndianCPU();

/**
 * Generic inplace endian order swapping.
 *
 * Be very careful, especially with "efficient" ones below
 * -> implementation library specific! Non-portable code!
 */

template<typename T>
T endianSwap(T& val);

/**
 * Some specific ones for efficiency
 */

float endianSwap(float& val);

int endianSwap(int& val);

unsigned endianSwap(unsigned& val);

short endianSwap(short& val);

} // namespace prlite

// Template include
#include "prlite_endian.tcc"

#endif // PRLITE_ENDIAN_HPP
