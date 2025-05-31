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

/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */


#ifndef PRLITE_GENDEFS_HPP
#define PRLITE_GENDEFS_HPP

// standard headers
#include <ciso646>  // and etc

 //namespace prlite{

// an alias for the floating point type used for model parameters, scores etc.
#ifndef REAL
# define REAL double // float //
#endif // REAL

// Slashes between directory names (Unix vs. Windows)
#ifdef _WIN32
#  define DIRSLASH '\\'
#endif // _WIN32


// Ensure __UNIX__ is defined if __linux__ is defined
#ifdef __linux__
# ifndef __UNIX__
#  define __UNIX__ __linux__
# endif
#endif

// Ensure __UNIX__ is defined if __APPLE__ is defined
#ifdef __APPLE__
# ifndef __UNIX__
#  define __UNIX__ __APPLE__
# endif
#endif

// under mingw you could define __UNIX__ to get / again.
#ifdef __UNIX__
#  undef DIRSLASH
#  define DIRSLASH '/'
#endif // __UNIX__

// __declspec(dllexport) exports all of the public data members and
// member functions in a class
#ifdef PRLITE_DLL
#ifdef PRLITE_EXPORTS
#define PRLITE_API __declspec(dllexport)
#else
#define PRLITE_API __declspec(dllimport)
#endif
#else
#define PRLITE_API
#endif

/*
 * The default setting used in constructors/configurers and reset members of models that can be
 * trained for the notUpdated status of model.
 *
 * false: paramUpdate will try to update the model whether it received training data or not.
 * Typically this causes a model that have not received any training data to bitch about this.
 *
 * true: typically it is not possible to distinguish between models that were already updated
 * and those that just received no training data. This makes everything a bit quieter. Since
 * higher level routines clean up redundant classes, the end result is pretty much the same.
 */
#define PRLITE_UPDATEDEF true

typedef unsigned char uchar;

//} // namespace prlite

#endif // PRLITE_GENDEFS_HPP
