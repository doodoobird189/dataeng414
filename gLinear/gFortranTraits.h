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

#ifndef GFORTRAN_TRAITS_H
#define GFORTRAN_TRAITS_H

// glinear headers
#include "gLinearMacros.h"

// standard headers
#include <complex>  // complex

/**@name Fortran Support Functions and types.
 * {\bf Note: } The following are all in the {\bf ::fortran} namespace.
 */
//@{
namespace fortran {
  /*
   * HaX0r'd straight from g2c.h shipped with egcs-2.91.66
   */
  /* F2C_INTEGER will normally be `int' but would be `long' on
   * 16-bit systems */
  /* we assume short, float are OK */
  /**
   * int
   */
  typedef int integer;
  /**
   * unsigned int
   */
  typedef unsigned int /* unsigned int */ uinteger;
  /**
   * address (pointer)
   */
  typedef char *address;
  /**
   * short int
   */
  typedef short int shortint;
  /**
   * floating point - single precision
   */
  typedef float real;
  /**
   * floating point - double precision
   */
  typedef double doublereal;
  /* more hacking by <goof@dsp.sun.ac.za>. this is to eliminate the name
   * clash between this complex and the one in the standard libs.
   */
  /**
   * complex floating point - single precision
   */
  struct complex {
    real r, i;
  }; // complex
  typedef struct complex complex;
  /**
   * complex floating point - double precision
   */
  typedef struct {
    doublereal r, i;
  } doublecomplex;
  /**
   * logical (bool) type
   */
  typedef long int /* long int */ logical;
  /**
   * logical (bool) type
   */
  typedef short int shortlogical;
  typedef char logical1;
  typedef char integer1;

  /* ANSI gets slightly upset about these two.
   */
  // typedef long long int /* long long */ longint;    /* system-dependent */
  // typedef unsigned long long int /* long long */ ulongint;  /* system-dependent */

# define TRUE_ (1)
# define FALSE_ (0)

  typedef long int /* long int */ flag;
  typedef long int /* long int */ ftnlen;
  typedef long int /* long int */ ftnint;


  /**
   * Object representing a Fortran array. This just provides a clean way
   * of generating arrays, with unknown length at compile-time, that we can
   * pass into Lapack routines. ANSI/ISO C++ does not like:
   *      type var[N];
   * decls if N is not a constant. These can then be replaced by:
   *     fortran_array<type> var(N);
   * which is not too messy.
   * The T* cast operator
   */
  template<typename T>
  class fortran_array {
  private:
    T* data;
    /**
     * Private default ctor. \\
     * It is impossible to call this. The size of this array is not supposed
     * to change and copy ctor and assignment is shallow, so it doesn't make
     * sense creating an empty one. \\
     */
    fortran_array(void);
  public:
    /**
     * Constructor specifying size. \\
     * As you probably would have gathered by now, this size can (and does
     * not) change after construction. \\
     * {\bf Warning: } No checking is done on the size parameter.
     * So please don't go do silly stuff like specifying a size of zero.
     * @param size Size of array (elements).
     */
    fortran_array(gLinear::gLinearTypes::T_Index size);

    /**
     * Copy constuctor. \\
     * Constructs a new object that is a {\bf shallow} copy of the argument
     * array. Shallow implies that it references the same storage. Please
     * read all the documentation of this class, especially the about the
     * destructor. It might save you some tears.
     * @param src Source array to copy.
     */
    fortran_array(const fortran_array& src);

    /**
     * Destructor. \\
     * Deallocate storage. The contents can be pretty much considered to
     * be gone after this call returns. \\
     * {\bf Warning: } Although it is possible to create more than one
     * {\em fortran_array} object pointing to the same storage, the first
     * one to destruct takes the storage with it (sorry, no reference
     * counting). This object should be treated like any other {\em T* v[N]}.
     */
    ~fortran_array(void);

    /**
     * Cast to T* operator. \\
     * Return pointer to the start of the data.
     * @return Start address of the array data as T*.
     */
    operator T*(void);

    /**
     * Cast to cosnt T* operator. \\
     * Return const pointer to the start of the data.
     * @return Start address of the array data as const T*.
     */
    operator const T*(void) const;

  }; // fortran_array
} // namespace fortran

GLINEAR_NAMESPACE_BEGIN

/**
 * Traits class for the Fortran types.
 */
template<typename T>
struct gFortranTraits {
  typedef void T_Real;
  typedef void T_Complex;
}; // gFortranTraits

/**
 * fortran::real
 */
template<>
struct gFortranTraits<fortran::real> {
  /**
   * Associated Fortran type
   */
  typedef fortran::real T_Fortran;
  /**
   * Asociated real type
   */
  typedef fortran::real T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::complex T_Complex;
  /**
   * Associated C++ type
   */
  typedef fortran::real T_Cpp;
}; // gFortranTraits
/**
 * fortran::complex
 */
template<>
struct gFortranTraits<fortran::complex> {
  /**
   * Associated Fortran type
   */
  typedef fortran::complex T_Fortran;
  /**
   * Asociated real type
   */
  typedef fortran::real T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::complex T_Complex;
  /**
   * Associated C++ type
   */
  typedef fortran::real T_Cpp;
}; // gFortranTraits
/**
 * fortran::doublereal
 */
template<>
struct gFortranTraits<fortran::doublereal> {

  /**
   * Associated Fortran type
   */
  typedef fortran::doublereal T_Fortran;
  /**
   * Asociated real type
   */
  typedef fortran::doublereal T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::doublecomplex T_Complex;
  /**
   * Associated C++ type
   */
  typedef std::complex<fortran::doublereal> T_Cpp;
}; // gFortranTraits
/**
 * fortran::doublecomplex
 */
template<>
struct gFortranTraits<fortran::doublecomplex> {
  /**
   * Associated Fortran type
   */
  typedef fortran::doublecomplex T_Fortran;
  /**
   * Asociated real type
   */
  typedef fortran::doublereal T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::doublecomplex T_Complex;
  /**
   * Associated C++ type
   */
  typedef std::complex<fortran::doublereal> T_Cpp;
}; // gFortranTraits

#if 0 // Apperently these are not needed (yet?) and they confuse MSC
/**
 * complex<fortran::real>
 */
template<>
struct gFortranTraits<std::complex<fortran::real> > {
  /**
   * Associated Fortran type
   */
  typedef fortran::complex T_Fortran;
  /**
   * Associated real type
   */
  typedef fortran::real T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::complex T_Complex;
}; // gFortranTraits
/*
 * complex<fortran::doublereal>
 */
template<>
struct gFortranTraits<std::complex<fortran::doublereal> > {
  /**
   * Associated fortran type
   */
  typedef fortran::doublecomplex T_Fortran;
  /**
   * Associated real type
   */
  typedef fortran::doublereal T_Real;
  /**
   * Associated complex type
   */
  typedef fortran::doublecomplex T_Complex;
}; // gFortranTraits
#endif // 0


//@}

GLINEAR_NAMESPACE_END

#include "gFortranTraits.tcc"

#endif // GFORTRAN_TRAITS_H
