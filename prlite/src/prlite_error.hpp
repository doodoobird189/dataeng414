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

/*******************************************************************************

   AUTHOR:  Edward de Villiers (Copyright University of Stellenbosch, all rights reserved.)
   DATE:    09/09/2003
   PURPOSE: A macro to report fatal errors and exit the program.  The parameter
            is an error message.

*******************************************************************************/

#ifndef PRLITE_ERROR_HPP
#define PRLITE_ERROR_HPP

// patrec headers

// standard headers
#include <iostream>  // cerr, endl
#include <sstream>  // ostringstream
#include <ciso646>  // and or not etc.

#ifdef ABORT_NOT_THROW
#include <signal.h>  // raise, SIGTRAP
#endif // ABORT_NOT_THROW

 //namespace prlite{

/**
 * Thrown when an ASSERT fails.
 *
 * ASSERTs check for conditions that may never happen.
 *
 * Therefore, if you feel you need to catch one of these, then the
 * code you're calling probably has to rather throw something else.
 * It is better to make that code use CHECK instead of ASSERT, and
 * rather throw a more appropriate error.
 */
class AssertFailed : public std::exception {
  std::string msg;
 public:
  AssertFailed(const std::string& _m = "") throw() : msg(_m) { }
  virtual ~AssertFailed() throw() { }
  virtual const char* what() const throw() { return msg.c_str(); }
};


/**
 * errors that higher level code may be expected to handle, or inform
 * a user about.  Do not throw this directly, rather see if there is a
 * derived exception class that fit the category of error your code
 * has found.
 *
 * Library specific exceptions should be defined in the appropriate
 * header file with the macro DEFINE_EXCEPTION() below.
 */
class PatRecError : public std::exception {
  std::string msg;
 public:
  PatRecError(const std::string& _m = "") throw() : msg(_m) { }
  PatRecError(const PatRecError& e) throw() : std::exception(e), msg(e.msg) { }
  virtual ~PatRecError() throw() { }
  virtual const char* what() const throw() { return msg.c_str(); }
};


/**
 * This is how we define the exception hierarchy.  PARENT should be
 * derived from PatRecError.
 */
#define PRLITE_DEFINE_EXCEPTION(TYPE, PARENT)                        \
  class TYPE : public PARENT {                                \
   public:                                                    \
   TYPE(const std::string& _m = "") throw() : PARENT(_m) { }  \
   TYPE(const PARENT& parent) throw() : PARENT(parent) { }    \
   virtual ~TYPE() throw() { }                                \
  }


/**
 * Thrown for IO errors: could not open a file, could not read from a
 * stream when more data is expected etc.
 */
PRLITE_DEFINE_EXCEPTION(IOError, PatRecError);

/**
 * Thrown for a special type of IO errors: file format violations.
 */
PRLITE_DEFINE_EXCEPTION(FileFormatError, IOError);

/**
 * Values out of bounds.  Mostly for checking things that could come
 * from user input.
 */
PRLITE_DEFINE_EXCEPTION(ValueError, PatRecError);


/*********************************************************
    Handy Macros
*/

/* === For fun ;-) ===
 *
 * Definitions of ANSI color codes:
 * <ESC>[<mode>;<color>m
 *
 * mode:
 *   00 = normal
 *   01 = bold
 * color:
 *   30 = black(normal) / dark gray(bold)
 *   31 = red
 *   32 = green
 *   33 = brown(normal) / yellow(bold)
 *   34 = blue
 *   35 = magenta
 *   36 = cyan
 *   37 = gray(normal) / white(bold)
 *
 */


/**
 * throw the given object as an exception.  this is preferred over using
 * "throw", since here it can be intercepted by ABORT_NOT_THROW.  This
 * makes getting a core dump easier.
 *
 * (remember ulimit -c unlimited!!!)
 */
#ifdef ABORT_NOT_THROW
#define PRLITE_THROW(E,M)                                                      \
  {                                                                     \
    std::cerr << '\n' << "[01;31mERROR =>[00m " << #E << ": " << M << std::endl; \
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;              \
    raise(SIGTRAP);                                                     \
  }
#else
#define PRLITE_THROW(E,M)                                                \
  {                                                               \
    std::ostringstream oss;                                       \
    oss << __FILE__ << ":" << __LINE__ << " " << #E << ": " << M; \
    throw E( oss.str() );                                         \
  }
#endif


/**
 * FAIL macro.  This is in the process of being deprecated.  Rather
 * use CHECK or ASSERT, since the semantics are better defined.
 */
#define PRLITE_FAIL(M) PRLITE_ASSERT(false, M)


/**
 * checks a condition C, else throws an exception of type E, with
 * message F.
 */
#define PRLITE_CHECK(C,E,F)                              \
  if ( ! (C) ) {                                  \
    PRLITE_THROW(E, "Check " << #C << " failed: " << F); \
  }


/**
 * asserts the condition C, else throws AssertFailed exception.
 *
 * Use ASSERT to check for conditions that may never ever happen,
 * because your code assumes or requires they don't.  If this fails,
 * the program should stop running.
 *
 * If you feel that some higher level piece of code that called your
 * code should be able to handle this, Rather call CHECK with an
 * appropriate exception type. (See the DEFINE_EXCEPTION calls above.)
 */
#define PRLITE_ASSERT(C,F) PRLITE_CHECK(C, AssertFailed, F)


//} // namespace prlite

#endif // PRLITE_ERROR_HPP

// This should be removed as when/if all the CAUTION calls have been
// changed to LOGWARN.  We don't really want other dependencies in
// prlite_error.hpp

// #include "prlite_logging.hpp"  // getPatrecLogger, Logger
