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

#ifndef GLINEAR_MACROS_H
#define GLINEAR_MACROS_H


// Static polymorphism
#define GLINEAR_STATIC_VIRTUAL(x) (static_cast<x*>(this))->

// DLL stuff for Windows: define GLINEAR_DLL, and GLINEAR_EXPORTS only
// when compiling gLinear DLL code.
#ifdef GLINEAR_DLL
#ifdef GLINEAR_EXPORTS
#define GLINEAR_API __declspec(dllexport)
#else
#define GLINEAR_API __declspec(dllimport)
#endif
#else
#define GLINEAR_API
#endif

// Namespace control
// this should be one of the following:
#define GLINEAR_USE_NAMESPACES // to put all the gLinear stuff in the gLinear namespace.
//  #undef GLINEAR_USE_NAMESPACES to not do it.
//  #undef GLINEAR_USE_NAMESPACES

#ifdef GLINEAR_USE_NAMESPACES
#  define GLINEAR_NAMESPACE_BEGIN namespace gLinear {
#  define GLINEAR_NAMESPACE_END }
#else
#  define GLINEAR_NAMESPACE_BEGIN
#  define GLINEAR_NAMESPACE_END
#endif

GLINEAR_NAMESPACE_BEGIN

/**
 * Tag class indicating that no initialization.
 */
struct gLinearNoInit {
}; // gLinearNoInit

/**
 * Const override macro to implement const policy.
 * @param x Type to which this should be cast.
 * @return De-const'ed this pointer.
 */
#define GLINEAR_CONST_OVERRIDE_THIS(x) const_cast<x*>(this)

/**
 * Class containing a number of global typedefs.
 */
class gLinearTypes {
public:
  /**
   * Type of index.
   */
  typedef int T_Index;
}; // gLinearTypes

/**
 * Tag classes
 */
struct GLINEAR_API gLinearTagBase {
}; // gLinearTagBase

/**
 * Tag class indicating a container implementing dense storage.
 */
struct GLINEAR_API DENSE : public gLinearTagBase {
  static const char className[];
}; // DENSE

/**
 * Tag class indicating a container implementing slice storage.
 */
struct GLINEAR_API SLICE : public gLinearTagBase {
  static const char className[];
}; // SLICE


/**
 * Tag class indicating a constainer implementing dense row-major storage.
 */
struct GLINEAR_API ROW_DENSE : public gLinearTagBase {
  static const char className[];
}; // ROW_DENSE

/**
 * Tag class indicating a constainer implementing dense column-major storage.
 */
struct GLINEAR_API COL_DENSE : public gLinearTagBase {
  static const char className[];
}; // COL_DENSE

/**
 * Tag class indicating a constainer implementing slice row-major storage.
 */
struct GLINEAR_API ROW_SLICE : public gLinearTagBase {
  static const char className[];
}; // ROW_SLICE

/**
 * Tag class indicating a constainer implementing slice column-major storage.
 */
struct GLINEAR_API COL_SLICE : public gLinearTagBase {
  static const char className[];
}; // COL_SLICE

/**
 * Tag class indicating a row.
 */
struct GLINEAR_API ROW : public gLinearTagBase {
  static const char className[];
}; // ROW

/**
 * Tag class indicating a column.
 */
struct GLINEAR_API COL : public gLinearTagBase {
  static const char className[];
}; // COL

GLINEAR_NAMESPACE_END

#endif // GLINEAR_MACROS_H
