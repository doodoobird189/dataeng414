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

// standard headers
#include <algorithm>  // transform
#include <complex>  // conj, real, imag, abs, arg, norm, cos, cosh, sin, sinh, exp, log, sqrt
#include <cmath>  // cos, cosh, sin, sinh, exp, log, sqrt


GLINEAR_NAMESPACE_BEGIN

template<typename T_Functor>
inline funcTransformVecInto<T_Functor>::funcTransformVecInto(const T_Functor& new_functor) :
  functor(new_functor) {
} // constructor

template<typename T_Functor>
template<typename T_Vector, typename T_Result>
inline void funcTransformVecInto<T_Functor>::operator() (const T_Vector& vec,
               T_Result result) {
  std::transform(vec.begin(), vec.end(), result.begin(), functor);
} // operator()

template<typename T_Functor>
template<typename T_Vector1, typename T_Vector2, typename T_Result>
inline void funcTransformVecInto<T_Functor>::operator() (const T_Vector1& vec1,
               const T_Vector2& vec2,
               T_Result result) {
  std::transform(vec1.begin(), vec1.end(), vec2.begin(), result.begin(), functor);
} // operator()

template<typename T_Value>
inline T_Value funcMinus<T_Value>::operator() (const T_Value& val) {
  return -val;
} // operator()

template<typename T_Value1, typename T_Value2>
inline funcAdd<T_Value1, T_Value2>::funcAdd(void)
  : constant(0) {
} // constructor

template<typename T_Value1, typename T_Value2>
inline funcAdd<T_Value1, T_Value2>::funcAdd(const T_Value2& new_constant) :
  constant(new_constant) {
} // constructor

template<typename T_Value1, typename T_Value2>
inline typename funcAdd<T_Value1, T_Value2>::T_Value
funcAdd<T_Value1, T_Value2>::operator() (const T_Value1& val) {
  return val + constant;
} // operator()

template<typename T_Value1, typename T_Value2>
inline typename funcAdd<T_Value1, T_Value2>::T_Value
funcAdd<T_Value1, T_Value2>::operator() (const T_Value1& val1,
           const T_Value2& val2) {
  return val1 + val2;
} // operator()

template<typename T_Value1, typename T_Value2>
inline typename funcSubtractCore<T_Value1, T_Value2, 0>::T_Value
funcSubtractCore<T_Value1, T_Value2, 0>::subtract(const T_Value1 first,
              const T_Value2 second) {
  return first - second;
} // subtract

template<typename T_Value1, typename T_Value2>
inline typename funcSubtractCore<T_Value1, T_Value2, 1>::T_Value
funcSubtractCore<T_Value1, T_Value2, 1>::subtract(const T_Value1 first,
              const T_Value2 second) {
  return second - first;
} // subtract

template<typename T_Value1, typename T_Value2, int Inv>
inline funcSubtract<T_Value1, T_Value2, Inv>::funcSubtract(void)
  : constant(0) {
} // constructor

template<typename T_Value1, typename T_Value2, int Inv>
inline funcSubtract<T_Value1, T_Value2, Inv>::funcSubtract(const T_Value2& new_constant)
  : constant(new_constant) {
} // constructor

template<typename T_Value1, typename T_Value2, int Inv>
inline typename funcSubtract<T_Value1, T_Value2, Inv>::T_Value
funcSubtract<T_Value1, T_Value2, Inv>::operator() (const T_Value1& val) {
  return funcSubtractCore<T_Value1, T_Value2, Inv>::subtract(val, constant);
} // operator()

template<typename T_Value1, typename T_Value2, int Inv>
inline typename funcSubtract<T_Value1, T_Value2, Inv>::T_Value
funcSubtract<T_Value1, T_Value2, Inv>::operator() (const T_Value1& val1,
               const T_Value2& val2) {
  return funcSubtractCore<T_Value1, T_Value2, Inv>::subtract(val1, val2);
} // operator()

template<typename T_Value1, typename T_Value2>
inline funcMultiply<T_Value1, T_Value2>::funcMultiply(void)
  : constant(1) {
} // constructor

template<typename T_Value1, typename T_Value2>
inline funcMultiply<T_Value1, T_Value2>::funcMultiply(const T_Value2& new_constant)
  : constant(new_constant) {
} // constructor

template<typename T_Value1, typename T_Value2>
inline typename funcMultiply<T_Value1, T_Value2>::T_Value
funcMultiply<T_Value1, T_Value2>::operator() (const T_Value1& val) {
  return val * constant;
} // operator()

template<typename T_Value1, typename T_Value2>
inline typename funcMultiply<T_Value1, T_Value2>::T_Value
funcMultiply<T_Value1, T_Value2>::operator() (const T_Value1& val1,
                const T_Value2& val2) {
  return val1 * val2;
} // operator()

template<typename T_Value1, typename T_Value2>
inline typename funcDivideCore<T_Value1, T_Value2, 0>::T_Value
funcDivideCore<T_Value1, T_Value2, 0>::divide(const T_Value1 first,
                const T_Value2 second) {
  return first / second;
} // divide

template<typename T_Value1, typename T_Value2>
inline typename funcDivideCore<T_Value1, T_Value2, 1>::T_Value
funcDivideCore<T_Value1, T_Value2, 1>::divide(const T_Value1 first,
                const T_Value2 second) {
  return second / first;
} // divide

template<typename T_Value1, typename T_Value2, int Inv>
inline funcDivide<T_Value1, T_Value2, Inv>::funcDivide(void)
  : constant(1) {
} // constructor

template<typename T_Value1, typename T_Value2, int Inv>
inline funcDivide<T_Value1, T_Value2, Inv>::funcDivide(const T_Value2& new_constant)
  : constant(new_constant) {
} // constructor

template<typename T_Value1, typename T_Value2, int Inv>
inline typename funcDivide<T_Value1, T_Value2, Inv>::T_Value
funcDivide<T_Value1, T_Value2, Inv>::operator() (const T_Value1& val) {
  return funcDivideCore<T_Value1, T_Value2, Inv>::divide(val, constant);
} // operator()

template<typename T_Value1, typename T_Value2, int Inv>
inline typename funcDivide<T_Value1, T_Value2, Inv>::T_Value
funcDivide<T_Value1, T_Value2, Inv>::operator() (const T_Value1& val1,
             const T_Value2& val2) {
  return funcDivideCore<T_Value1, T_Value2, Inv>::divide(val1, val2);
} // operator()

template<typename T_Vector1, typename T_Vector2>
inline typename funcDot<T_Vector1, T_Vector2>::T_Value
funcDot<T_Vector1, T_Vector2>::dodot(const T_Vector1& vec1,
             const T_Vector2& vec2) {
  T_Value result = 0;
  const typename T_Vector2::T_ConstIterator iter2_end = vec2.end();
  typename T_Vector1::T_ConstIterator iter1 = vec1.begin();
  for (typename T_Vector2::T_ConstIterator iter2 = vec2.begin();
       iter2 != iter2_end;
       ++iter2, ++iter1) {
    result += (*iter2) * (*iter1);
  } // for
  return result;
} // dodot

template<typename T_Vector1, typename T_Vector2>
inline funcDot<T_Vector1, T_Vector2>::funcDot(void) {
} // default constructor

template<typename T_Vector1, typename T_Vector2>
inline funcDot<T_Vector1, T_Vector2>::funcDot(const T_Vector2& new_constant)
  : constant(new_constant) {
} // constructor

template<typename T_Vector1, typename T_Vector2>
inline typename funcDot<T_Vector1, T_Vector2>::T_Value funcDot<T_Vector1, T_Vector2>::operator() (const T_Vector1& vec1) {
  return dodot(vec1, constant);
} // operator()

template<typename T_Vector1, typename T_Vector2>
inline typename funcDot<T_Vector1, T_Vector2>::T_Value
funcDot<T_Vector1, T_Vector2>::operator() (const T_Vector1& vec1,
             const T_Vector2& vec2) {
  return dodot(vec1, vec2);
} // operator()

template<typename T>
inline typename funcConj<T>::T_Complex
funcConj<T>::operator() (const T_Complex& val) {
  return std::conj(val);
} // operator()

template<typename T>
inline typename funcConj<T>::T_Real
funcConj<T>::operator() (const T_Real& val) {
  return val;
} // operator()

template<typename T>
inline typename funcReal<T>::T_Real
funcReal<T>::operator() (const T_Complex& val) {
  return std::real(val);
} // operator()

template<typename T>
inline typename funcReal<T>::T_Real
funcReal<T>::operator() (const T_Real& val) {
  return val;
} // operator()

template<typename T>
inline typename funcImag<T>::T_Real
funcImag<T>::operator() (const T_Complex& val) {
  return std::imag(val);
} // operator()

template<typename T>
inline typename funcImag<T>::T_Real
funcImag<T>::operator() (const T_Real& val) {
  return 0;
} // operator()

template<typename T>
inline typename funcAbs<T>::T_Real
funcAbs<T>::operator() (const T& val) {
  return std::abs(val);
} // operator()

template<typename T>
inline typename funcArg<T>::T_Real funcArg<T>::operator() (const T_Complex& val) {
  return std::arg(val);
} // operator()

template<typename T>
inline typename funcArg<T>::T_Real funcArg<T>::operator() (const T_Real& val) {
  return 0;
} // operator()

template<typename T>
inline typename funcNorm<T>::T_Real funcNorm<T>::operator() (const T_Complex& val) {
  return std::norm(val);
} // operator()

template<typename T>
inline typename funcNorm<T>::T_Real funcNorm<T>::operator() (const T_Real& val) {
  return val * val;
} // operator()

template<typename T>
inline T funcCos<T>::operator() (const T& val) {
  return std::cos(val);
} // operator()

template<typename T>
inline T funcCosh<T>::operator() (const T& val) {
  return std::cosh(val);
} // operator()

template<typename T>
inline T funcSin<T>::operator() (const T& val) {
  return std::sin(val);
} // operator()

template<typename T>
inline T funcSinh<T>::operator() (const T& val) {
  return std::sinh(val);
} // operator()

template<typename T>
inline T funcExp<T>::operator() (const T& val) {
  return std::exp(val);
} // operator()

template<typename T>
inline T funcLog<T>::operator() (const T& val) {
  return std::log(val);
} // operator()

template<typename T>
inline T funcSqrt<T>::operator() (const T& val) {
  return std::sqrt(val);
} // operator()

template<typename T>
inline T funcStep<T>::operator() (const T& val) {
  return (val > static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(0);
} // operator()

template<typename T>
inline T funcImpulse<T>::operator() (const T& val) {
  return (val == static_cast<T>(0)) ? static_cast<T>(1) : static_cast<T>(0);
} // operator()

GLINEAR_NAMESPACE_END
