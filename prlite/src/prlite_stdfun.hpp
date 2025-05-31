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
 * Author        | Dr J.A. du Preez, et. al.                                  *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 1995/09/14                                                 *
 * Description   | Various useful routines                                    *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

#ifndef PRLITE_STDFUN_HPP
#define PRLITE_STDFUN_HPP

// patrec headers
#include "prlite_genvec.hpp"  // Vector
#include "prlite_gendefs.hpp"  // REAL

// standard headers
#include <math.h>

/**************************** Constants & Defines      ***********************/

namespace prlite{

PRLITE_API extern const float INFTY;

#ifndef M_PI
#  define M_PI 3.141592653589793115997963468544185161590576171875
#endif

/************************** Type conversion (rounding) ***********************/

/**
 * Round floating-point number to nearest integer.
 * Overloaded for various floating-point and integer variants.
 */
template<typename T>
int round2int(T x);

template<typename T>
long round2long(T x);

template<typename T>
unsigned round2unsigned(T x);

template<typename T>
unsigned long round2ulong(T x);

/**
 * Rounds to a specified precision of decimal points.
 * E.g. fround(31.1415, 1) -> 31.1
 * Also works with precision < 0 eg fround(31.1415, -1) -> 30.0
 */
extern float fround(float x,
                    int prec);
extern int truncFloatToInt(float x);

/************************** Math functions             ***********************/

extern float cln(float);
// extern Vector<float> log(const Vector<float>& arr);
/// return log(e^L1 + e^L2)

extern int sign(float x);
extern float absfloat(float);
extern double absdouble(double);
extern int even(long);

// also see almostEqual in prlite_testing.hpp
template<typename T>
inline bool realEqual(T x,
                      T y);

/**
 * Projects x on the line connecting (x1, y1) and (x2, y2).
 */
extern float linInterpolate(float x1,
                            float y1,
                            float x2,
                            float y2,
                            float x);

/**
 * Least-squares polynomial curve fitting.
 * (RvdM 18/09/97)
 * @param x independent variable sequence
 * @param y dependent variable sequence
 * @param m order of polynomial
 * @param coef polynomial coefficients [a_0 , a_1 , ... , a_m] (return value)
 * @return the polynomial regression error
 */
extern float polyFit(const gLinear::gRowVector<float>& x,
                     const gLinear::gRowVector<float>& y,
                     int m,
                     gLinear::gRowVector<float>& coef);

/**
 * Polynomial evaluation.
 * Evaluate a polynomial specified by the coefficients
 * [a_0, a_1, ... , a_m] in coef at the values of x
 */
extern gLinear::gRowVector<float> polyEval(const gLinear::gRowVector<float>& x,
            const gLinear::gRowVector<float>& coef);

/**
 * Does a rapid search in a Vector to find the position of a specific
 * input value, then returns the corresponding lookup/output
 * value. Can interpolate, but does not extrapolate (ie when the input
 * goes outside the given input range the return value will be fixed
 * at the current outer limits). Useful for amongst others, PAV
 * lookup.
 *
 * @param x The input value for which a corresponding lookup/output
 * value is required.
 *
 * @param sortedX The series of known lookup input values. Should be
 * strictly increasing with no equal values.
 *
 * @param outVals The corresponding lookup values for each of the
 * sortedX values. Must therefore have same length as sortedX.
 *
 * @param interp Governs behaviour if x falls between two successive
 * values in sortedX. If 0 the outVals corresponding to the smaller
 * sortedX will be returned (ie a left-continuous step curve is
 * assumed for outVals). If not 0 the result will be a linear
 * interpolation of the two possible output values. NOT DONE: Should
 * be easy to extend to higher order interpolations in which case
 * interp will give the order of the interpolation
 * polinomial. Currently any value > 0 will give linear interpolation.
 */
extern float sortedLookup( float x,
                           const gLinear::gRowVector<float>& sortedX,
                           const gLinear::gRowVector<float>& outVals,
                           unsigned interp);


/**
 * Calculate "a combination b".
 */
template<typename T>
T combination(int a, int b);

/**
 * Calculates the next power of 2 equal or above to 'from'.
 */
extern int nextPowerOf2(unsigned from);

/**
 * Returns the non-negative integer exponent of an integer.
 */
extern int intPow(const int base,
                  const unsigned int exponent);

/************************** Frame indexing             ***********************/

/**
 * Returns the frameNo at a specific end-time i.e. the frameNo i
 * for which the time t is in the range ( t_mid(i) - Ts/2 .. t_mid(i) + Ts/2].
 * Remember that this is endtimes, at t = 0 you will get frame -1 or earlier.
 */
extern int frameIdx( double t,
                     double frameShift );

/**
 * Returns the end-time of frame defined as halfway between the
 * midpoint of this frame and the midpoint of the next.
 */
extern double frameTime( int idx,
                         double frameShift );

extern unsigned noOfVectorsInFrames(unsigned noOfFrames,
            unsigned frmLen,
            unsigned frmShift);

} // namespace prlite

#include "prlite_stdfun.tcc"

#endif // PRLITE_STDFUN_HPP
