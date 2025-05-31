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

#ifndef GVECTOR_OPERS_H
#define GVECTOR_OPERS_H

// glinear headers
#include "gLinearMacros.h"
#include "gRowVector.h"
#include "gColVector.h"
#include "gLinearFunctors.h"
#include "gCppTraits.h"

// standard headers
#include <iosfwd>  // ostream, istream


GLINEAR_NAMESPACE_BEGIN


/**@name Unary Vector operators
 */
//@{
/**
 * gRowVector<T, T_Stor> unary minus.
 * Return vector with operator-(T) applied to each element.
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
operator-(const gRowVector<T, T_Stor>& vec);

/**
 * gColVector<T, T_Stor> unary minus.
 * Return vector with operator-(T) applied to each element.
 */
template<typename T, typename T_Stor>
typename gColVector<T, T_Stor>::T_Dense
operator-(const gColVector<T, T_Stor>& vec);

// support for operator-()
template<typename T_Vec>
typename T_Vec::T_Dense o__vecUnaryMinus(const T_Vec& vec);
//@}

/**@name gRowVector<T1, stortype> - gRowVector<T2, stortype> operators
 */
//@{
/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> equality
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator==(const gRowVector<T1, lhstortype>& lharg,
    const gRowVector<T2, rhstortype>& rharg);
/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> inequality
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator!=(const gRowVector<T1, lhstortype>& lharg,
    const gRowVector<T2, rhstortype>& rharg);

/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> lessthan.
 * This operator allows an ordering of gRowVectors, enabling their use
 *   in various other classes.
 * @return false if lharg.size() > rharg.size(),
 *  true if either lharg.size() < rharg.size() or
 *  if lharg.size() == rharg.size() an element by element
 *  comparison is done starting at index 0 to find the first
 *  varying pair which then determines the outcome.
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator<(const gRowVector<T1, lhstortype>& lharg,
         const gRowVector<T2, rhstortype>& rharg);

/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> addition
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gRowVector<T1, lhstortype>& lharg,
    const gRowVector<T2, rhstortype>& rharg);
/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, lhstortype> subtraction
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gRowVector<T1, lhstortype>& lharg,
    const gRowVector<T2, rhstortype>& rharg);
//@}

/**@name gColVector<T1, stortype> - gColVector<T2, stortype> operators
 */
//@{
/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> equality
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator==(const gColVector<T1, lhstortype>& lharg,
    const gColVector<T2, rhstortype>& rharg);
/**
 * gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> inequality
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator!=(const gColVector<T1, lhstortype>& lharg,
    const gColVector<T2, rhstortype>& rharg);

/**
 * gColVector<T1, lhstortype> - gColVector<T2, rhstortype> lessthan.
 * This operator allows an ordering of gColVectors, enabling their use
 *   in various other classes.
 * @return false if lharg.size() > rharg.size(),
 *  true if either lharg.size() < rharg.size() or
 *  if lharg.size() == rharg.size() an element by element
 *  comparison is done starting at index 0 to find the first
 *  varying pair which then determines the outcome.
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool operator<(const gColVector<T1, lhstortype>& lharg,
         const gColVector<T2, rhstortype>& rharg);

/**
 * gColVector<T1, lhstortype> - gColVector<T2, rhstortype> addition
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gColVector<T1, lhstortype>& lharg,
    const gColVector<T2, rhstortype>& rharg);
/**
 * gColVector<T1, lhstortype> - gColVector<T2, rhstortype> subtraction
 */
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gColVector<T1, lhstortype>& lharg,
    const gColVector<T2, rhstortype>& rharg);
//@}

/**@name gRowVector<T, stortype> - gColVector<T, stortype> operators
 */

//@{
/**
 * gRowVector<T, LHStor> - gColVector<T, RHStor> multiplication (dot product)
 */
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gNumPrec2<T1, T2>::T_Highest operator*(const gRowVector<T1, T_LHStor>& vecl,
            const gColVector<T2, T_RHStor>& vecr);

//@}


/**@name gRowVector<T, stortype> - Scalar
 */
//@{

//----------------------jadp start of mods  13/7/2001

/**
 * Vector == Scalar :  equality
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator==(const gRowVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector != Scalar :  inequality
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator!=(const gRowVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector > Scalar :  bigger
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator>(const gRowVector<T1, T_Stor>& lharg,
          const T2& rharg);

/**
 * Vector < Scalar :  smaller
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator<(const gRowVector<T1, T_Stor>& lharg,
          const T2& rharg);

/**
 * Vector >= Scalar :  equal or bigger
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator>=(const gRowVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector <= Scalar :  equal or smaller
 */
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE> operator<=(const gRowVector<T1, T_Stor>& lharg,
           const T2& rharg);

//----------------------jadp end of mods  13/7/2001

/**
 * Vector + Scalar :  addition
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gRowVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector - Scalar :  subtraction
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gRowVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector * Scalar :  multiplication
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const gRowVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector / Scalar :  division
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const gRowVector<T1, vecstortype>& vec,
    const T2& scalar);
//@}

/**@name Scalar - gRowVector<T, stortype>
 */
//@{
/**
 * Scalar - Vector addition
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const T1& scalar,
    const gRowVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector subtraction
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const T1& scalar,
    const gRowVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector multiplication
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const T1& scalar,
    const gRowVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector division
 */
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const T1& scalar,
    const gRowVector<T2, vecstortype>& vec);
//@}

/**@name gColVector<T, stortype> - Scalar
 */
//@{

//----------------------jadp start of mods  13/7/2001

/**
 * Vector == Scalar :  equality
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator==(const gColVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector != Scalar :  inequality
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator!=(const gColVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector > Scalar :  bigger
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator>(const gColVector<T1, T_Stor>& lharg,
          const T2& rharg);

/**
 * Vector < Scalar :  smaller
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator<(const gColVector<T1, T_Stor>& lharg,
          const T2& rharg);

/**
 * Vector >= Scalar :  bigger or equal
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator>=(const gColVector<T1, T_Stor>& lharg,
           const T2& rharg);

/**
 * Vector <= Scalar :  smaller or equal
 */
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE> operator<=(const gColVector<T1, T_Stor>& lharg,
           const T2& rharg);

//----------------------jadp end of mods  13/7/2001

/**
 * Vector + :  Scalar addition
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gColVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector - :  Scalar subtraction
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gColVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector * Scalar :  multiplication
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const gColVector<T1, vecstortype>& vec,
    const T2& scalar);
/**
 * Vector / Scalar :  division
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const gColVector<T1, vecstortype>& vec,
    const T2& scalar);
//@}

/**@name Scalar - gColVector<T, stortype>
 */
//@{
/**
 * Scalar - Vector addition
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const T1& scalar,
    const gColVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector subtraction
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const T1& scalar,
    const gColVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector multiplication
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const T1& scalar,
    const gColVector<T2, vecstortype>& vec);
/**
 * Scalar - Vector division
 */
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const T1& scalar,
    const gColVector<T2, vecstortype>& vec);
//@}

//----------------------jadp start of mods  13/7/2001

/**@name gRowVector<bool, T_Stor1> - gRowVector<bool, T_Stor2> operators
 */
//@{

/**
 * gRowVector<bool, lhstortype> && gRowVector<bool, rhstortype> :  logical and
 */
template<typename lhstortype, typename rhstortype>
gRowVector<bool, DENSE> operator&&(const gRowVector<bool, lhstortype>& lharg,
           const gRowVector<bool, rhstortype>& rharg);

/**
 * gRowVector<bool, lhstortype> || gRowVector<bool, rhstortype> :  logical or
 */
template<typename lhstortype, typename rhstortype>
gRowVector<bool, DENSE> operator||(const gRowVector<bool, lhstortype>& lharg,
           const gRowVector<bool, rhstortype>& rharg);

//@}

/**@name gColVector<bool, T_Stor1> - gColVector<bool, T_Stor2> operators
 */
//@{

/**
 * gColVector<bool, lhstortype> && gColVector<bool, rhstortype> :  logical and
 */
template<typename lhstortype, typename rhstortype>
gColVector<bool, DENSE> operator&&(const gColVector<bool, lhstortype>& lharg,
           const gColVector<bool, rhstortype>& rharg);

/**
 * gColVector<bool, lhstortype> || gColVector<bool, rhstortype> :  logical or
 */
template<typename lhstortype, typename rhstortype>
gColVector<bool, DENSE> operator||(const gColVector<bool, lhstortype>& lharg,
           const gColVector<bool, rhstortype>& rharg);

//@}

//----------------------jadp end of mods  13/7/2001

//----------------------lude start of mods  22/11/2001

/**@name Unary Vector Functions
 */
//@{

/**
 * Basic unary vector function template. \\
 * This differs from funcTransformVecInto by being a function instead of a class,
 * since that's the way functions typically get called. Although funcTransformVecInto
 * (and transform) is all nice and useful as a function with its operator(), it still
 * requires the user to instantiate the appropriate functor class(es) beforehand.
 * This template function removes the schlepp for those basic callables.
 * The output container is always dense.
 */
template <typename T, typename T_Stor, typename T_Functor>
typename gRowVector<T, T_Stor>::T_Dense
basicFunc(const gRowVector<T, T_Stor>& arg);

/**
 * Basic unary vector function template for complex->real. \\
 * This handles the case where a complex argument is reduced to a real number
 * (as in abs()). The output container is again dense. Maybe the DENSE spec could
 * be more generic/traited?
 */
template <typename T, typename T_Stor, typename T_Functor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
complexToRealFunc(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the complex conugate of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
conj(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the real components of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
real(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the imaginary components of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
imag(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the absolute value of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
abs(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the argument of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
arg(const gRowVector<T, T_Stor>& arg);

/**
 * Return vector containing the norm of the
 * argument vector elements
 */
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
norm(const gRowVector<T, T_Stor>& arg);

/**
 * Cosine
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
cos(const gRowVector<T, T_Stor>& arg);

/**
 * Hyperbolic Cosine
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
cosh(const gRowVector<T, T_Stor>& arg);

/**
 * Sine
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sin(const gRowVector<T, T_Stor>& arg);

/**
 * Hyperbolic Sine
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sinh(const gRowVector<T, T_Stor>& arg);

/**
 * Exponent
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
exp(const gRowVector<T, T_Stor>& arg);

/**
 * Natural logarithm
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
log(const gRowVector<T, T_Stor>& arg);

/**
 * Square root
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sqrt(const gRowVector<T, T_Stor>& arg);

/**
 * Unit step function
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
step(const gRowVector<T, T_Stor>& arg);

/**
 * Impulse (dirac-delta) function
 */
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
impulse(const gRowVector<T, T_Stor>& arg);

//@}

//----------------------lude end of mods  22/11/2001

/**@name Human-readable Stream I/O
 */
//@{
/**
 * Write vector length, followed by contents, to ostream
 */
template<typename T, typename vecstortype>
std::ostream& operator<<(std::ostream& os, const gRowVector<T, vecstortype>& vec);

/**
 * Write vector length, followed by contents, to ostream
 */
template<typename T, typename vecstortype>
std::ostream& operator<<(std::ostream& os, const gColVector<T, vecstortype>& vec);

/**
 * Read vector length, followed by contents, from istream
 */
template<class T, typename vecstortype>
std::istream& operator>>(std::istream& is, gRowVector<T, vecstortype>& vec);

/**
 * Read vector length, followed by contents, from istream
 */
template<typename T, typename vecstortype>
std::istream& operator>>(std::istream& is, gColVector<T, vecstortype>& vec);

//@}

/**@name Binary Stream I/O
 */
//@{
/**
 * Write vector length, followed by contents, to ostream (raw binary)
 */
template<typename T, typename vecstortype>
std::ostream& global_write(std::ostream& os, const gRowVector<T, vecstortype>& vec);

/**
 * Write vector length, followed by contents, to ostream (raw binary)
 */
template<typename T, typename vecstortype>
std::ostream& global_write(std::ostream& os, const gColVector<T, vecstortype>& vec);

/**
 * Read vector length, followed by contents, from istream (raw binary)
 */
template<class T, typename vecstortype>
std::istream& global_read(std::istream& is, gRowVector<T, vecstortype>& vec);

/**
 * Read vector length, followed by contents, from istream (raw binary)
 */
template<typename T, typename vecstortype>
std::istream& global_read(std::istream& is, gColVector<T, vecstortype>& vec);

//@}

//----------------------jadp start of mods  13/7/2001

/**@name Functions
 */
//@{

/**
 * Forms a (possibly shorter) vector of the indices where gRowVector vec is true
 */
template<typename T_Stor>
gRowVector<gLinearTypes::T_Index, DENSE> find(const gRowVector<bool, T_Stor>& vec);

/**
 * Forms a (possibly shorter) vector of the indices where gColVector vec is true
 */
template<typename T_Stor>
gColVector<gLinearTypes::T_Index, DENSE> find(const gColVector<bool, T_Stor>& vec);

//@}

//----------------------jadp end of mods  13/7/2001

GLINEAR_NAMESPACE_END

#include "gVectorOpers.tcc"

#endif // GVECTOR_OPERS_H
