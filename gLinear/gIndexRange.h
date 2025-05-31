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

#ifndef GINDEX_RANGE_H
#define GINDEX_RANGE_H

// glinear headers
#include "gLinearMacros.h"

GLINEAR_NAMESPACE_BEGIN

/**
 * Class storing subrange index subrange parameters.
 */
class GLINEAR_API gIndexRange {
public:
  typedef gLinearTypes::T_Index T_Index;
  /**
   * Index of the start of the subrange. \\
   * Zero or positive integer.
   */
  T_Index start;
  /**
   * Index of the end of the subrange. \\
   * Zero or positive integer.
   */
  T_Index end;
  /**
   * Stride value. \\
   * Zero or positive integer.
   */
  T_Index stride;
  /**
   * Construct single-valued range. \\
   * The start and end values are both set equal to the supplied argument
   * {\em val}.
   * The stride is set equal to 1.
   * @param val start and end value.
   */
  gIndexRange(T_Index val);
  /**
   * Construct proper subrange.
   * The start, end and stride values are initialized from the arguments.
   * @param n_start New start value.
   * @param n_end New end value.
   * @param n_stride New stride value.
   */
  gIndexRange(T_Index new_start,
              T_Index new_end,
              T_Index new_stride = 1);
  /**
   * Set subrange values.
   * The start, end and stride values are set from the arguments.
   * @param n_start New start value.
   * @param n_end New end value.
   * @param n_stride New stride value.
   */
  void setRange(T_Index new_start,
                T_Index new_end,
                T_Index new_stride);
  /**
   * Return start index value.
   */
  T_Index getStart(void) const;
  /**
   * Return end index value.
   */
  T_Index getEnd(void) const;
  /**
   * Return stride value.
   */
  T_Index getStride(void) const;

}; // gIndexRange

GLINEAR_NAMESPACE_END

#endif // GINDEX_RANGE_H
