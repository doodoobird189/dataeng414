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

#ifndef GLINEAR_FUNCTORS_H
#define GLINEAR_FUNCTORS_H

// glinear headers
#include "gCppTraits.h"

GLINEAR_NAMESPACE_BEGIN

/**
 * Base class for all functors
 */
struct gLinearFunctorBase {
}; // gLinearFunctorBase


/**
 * Apply functor to each vector element or pair of elements,
 * storing the result in {\em result}.
 */
template<typename T_Functor>
class funcTransformVecInto : public gLinearFunctorBase {
private:
  T_Functor functor;
public:
  funcTransformVecInto(const T_Functor& new_functor);

  // apply functor to each vector element, updating result
  template<typename T_Vector, typename T_Result>
  void operator() (const T_Vector& vec, T_Result result);

  // apply functor to each pair of vector elements, updating result
  template<typename T_Vector1, typename T_Vector2, typename T_Result>
  void operator() (const T_Vector1& vec1, const T_Vector2& vec2, T_Result result);

}; // funcTransformVecInto

/**
 * Return -val
 */
template<typename T_Value>
class funcMinus : public gLinearFunctorBase {
public:
  T_Value operator() (const T_Value& val);
}; // funcMinus

/**
 * Return val1 + val2 or val + constant.
 */
template<typename T_Value1, typename T_Value2 = T_Value1>
class funcAdd : public gLinearFunctorBase {
private:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  T_Value2 constant;
public:
  funcAdd(void);

  funcAdd(const T_Value2& new_constant);

  T_Value operator() (const T_Value1& val);

  T_Value operator() (const T_Value1& val1,
                      const T_Value2& val2);

}; // funcAdd


// primary decl.
template <typename T_Value1, typename T_Value2, int Inv>
class funcSubtractCore;

// non-inverting specialization
template<typename T_Value1, typename T_Value2>
class funcSubtractCore<T_Value1, T_Value2, 0> {
public:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  static T_Value subtract(const T_Value1 first, const T_Value2 second);
}; // funcSubtractCore

// inverting specialization
template<typename T_Value1, typename T_Value2>
class funcSubtractCore<T_Value1, T_Value2, 1> {
public:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  static T_Value subtract(const T_Value1 first, const T_Value2 second);
}; // funcSubtractCore

/**
 * Return val1 - val2 or val - constant.
 */
template<typename T_Value1, typename T_Value2 = T_Value1, int Inv = 0>
class funcSubtract : public gLinearFunctorBase {
private:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  T_Value2 constant;
public:
  funcSubtract(void);

  funcSubtract(const T_Value2& new_constant);

  T_Value operator() (const T_Value1& val);

  T_Value operator() (const T_Value1& val1,
                      const T_Value2& val2);

}; // funcSubtract


/**
 * Return val1 * val2 of val * constant.
 */
template<typename T_Value1, typename T_Value2 = T_Value1>
class funcMultiply : public gLinearFunctorBase {
private:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  T_Value2 constant;
public:
  funcMultiply(void);

  funcMultiply(const T_Value2& new_constant);

  T_Value operator() (const T_Value1& val);

  T_Value operator() (const T_Value1& val1,
                      const T_Value2& val2);

}; // funcMultiply


// primary decl.
template <typename T_Value1, typename T_Value2, int Inv>
class funcDivideCore;

// non-inverting specialization
template<typename T_Value1, typename T_Value2>
class funcDivideCore<T_Value1, T_Value2, 0> {
public:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  static T_Value divide(const T_Value1 first, const T_Value2 second);
}; // funcDivideCore

// inverting specialization
template<typename T_Value1, typename T_Value2>
class funcDivideCore<T_Value1, T_Value2, 1> {
public:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  static T_Value divide(const T_Value1 first, const T_Value2 second);
}; // funcDivideCore

/**
 * Return val1 / val2 or val / constant.
 */
template<typename T_Value1, typename T_Value2 = T_Value1, int Inv = 0>
class funcDivide : public gLinearFunctorBase {
private:
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
  T_Value2 constant;
public:
  funcDivide(void);

  funcDivide(const T_Value2& new_constant);

  T_Value operator() (const T_Value1& val);

  T_Value operator() (const T_Value1& val1,
                      const T_Value2& val2);

}; // funcDivide

/**
 * Return inner(vec1, vec2) or inner(vec, constant)
 */
template<typename T_Vector1, typename T_Vector2>
class funcDot : public gLinearFunctorBase {
  typedef typename T_Vector1::T_Value T_Value1;
  typedef typename T_Vector2::T_Value T_Value2;
  typedef typename gNumPrec2<T_Value1, T_Value2>::T_Highest T_Value;
private:
  T_Vector2 constant;
  T_Value dodot(const T_Vector1& vec1, const T_Vector2& vec2);

public:
  typedef T_Value T_Result;
  funcDot(void);

  funcDot(const T_Vector2& new_constant);

  T_Value operator() (const T_Vector1& vec1);

  T_Value operator() (const T_Vector1& vec1,
                      const T_Vector2& vec2);

}; // funcDot

/**
 * Return the complex conjugate
 */
template<typename T>
class funcConj : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Complex operator() (const T_Complex& val);

  T_Real operator() (const T_Real& val);
}; // funcConj

/**
 * Return the real component
 */
template<typename T>
class funcReal : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Real operator() (const T_Complex& val);

  T_Real operator() (const T_Real& val);

}; // funcReal

/**
 * Return the imaginary component
 */
template<typename T>
class funcImag : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Real operator() (const T_Complex& val);

  T_Real operator() (const T_Real& val);

}; // funcImag

/**
 * Return the absolute value
 */
template<typename T>
class funcAbs : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Real operator() (const T& val);

}; // funcAbs

/**
 * Return the argument (angle or arc tangent of imag / real)
 */
template<typename T>
class funcArg : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Real operator() (const T_Complex& val);

  T_Real operator() (const T_Real& val);

}; // funcArg

/**
 * Return the norm (real^2 + imag^2)
 */
template<typename T>
class funcNorm : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T_Real operator() (const T_Complex& val);

  T_Real operator() (const T_Real& val);

}; // funcNorm

/**
 * Return the cosine
 */
template<typename T>
class funcCos : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);

}; // funcCos

/**
 * Return the hyperbolic cosine
 */
template<typename T>
class funcCosh : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcCosh

/**
 * Return the sine
 */
template<typename T>
class funcSin : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcSin

/**
 * Return the hyperbolic sine
 */
template<typename T>
class funcSinh : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcSinh

/**
 * Return the exponent
 */
template<typename T>
class funcExp : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcExp

/**
 * Return the natural logarithm
 */
template<typename T>
class funcLog : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcLog

/**
 * Return the square root
 */
template<typename T>
class funcSqrt : public gLinearFunctorBase {
  typedef typename gCppTraits<T>::T_Complex T_Complex;
  typedef typename gCppTraits<T>::T_Real T_Real;
public:
  T operator() (const T& val);
}; // funcSqrt

/**
 * Unit step function: u(n)
 */
template<typename T>
class funcStep : public gLinearFunctorBase {
public:
  T operator() (const T& val);
}; // funcStep

/**
 * Impulse (dirac-delta) function
 */
template<typename T>
class funcImpulse : public gLinearFunctorBase {
public:
  T operator() (const T& val);
}; // funcImpulse

GLINEAR_NAMESPACE_END

#include "gLinearFunctors.tcc"

#endif // GLINEAR_FUNCTORS_H
