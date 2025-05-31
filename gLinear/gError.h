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


#ifndef GERROR_H
#define GERROR_H

// glinear headers
#include "gLinearMacros.h"

// standard headers
#include <string>  // string

GLINEAR_NAMESPACE_BEGIN

typedef gLinearTypes::T_Index T_Index;


/**
 * Error handler base class. \\
 * All error objects are derived from this class. These objects are
 * constructed when an error is encountered. The contructor then calls
 * a member function in this base class that handles the error.
 */
class GLINEAR_API gErrorBaseHandler {
protected:
  void handle(const std::string& errstr,
              int exit_code = -1);
}; // gErrorBaseHandler

/**
 * 1D and 2D index out of bounds
 */
class GLINEAR_API gErrorIndex : public gErrorBaseHandler {
public:
  /*
   * General 2D Index out of bounds
   */
  gErrorIndex(std::string classname,
              std::string funcname,
              T_Index maxrows,
              T_Index maxcols,
              T_Index row,
              T_Index col);
  /*
   * General 1D Index out of bounds
   */
  gErrorIndex(std::string classname,
              std::string funcname,
              std::string errname,
              T_Index maxval,
              T_Index val);
}; // gErrorIndex

/**
 * Slice error
 */
class GLINEAR_API gErrorSlice : public gErrorBaseHandler {
public:
  gErrorSlice(std::string classname,
              std::string funcname,
              std::string errname,
              T_Index maxindx,
              T_Index strt,
              T_Index end,
              T_Index stride);
}; // gErrorSlice

/**
 * 1D and 2D subrange out of bounds
 */
class GLINEAR_API gErrorSubrange : public gErrorBaseHandler {
public:
  /**
   * General 1D Subrange out of bounds
   */
  gErrorSubrange(std::string classname,
                 std::string funcname,
                 std::string errname,
                 T_Index from,
                 T_Index to,
                 T_Index maxval);
}; // gErrorSubrange

/**
 * Aliases exist
 */
class GLINEAR_API gErrorAliases : public gErrorBaseHandler {
public:
  gErrorAliases(std::string classname,
                std::string funcname,
                std::string errname);
}; // gErrorAliases

/**
 * Shape error
 */
class GLINEAR_API gErrorShape : public gErrorBaseHandler {
public:
  gErrorShape(std::string classname,
              std::string funcname,
              std::string errname,
              T_Index rows,
              T_Index cols);
}; // gErrorShape

/**
 * Size mismatch
 */
class GLINEAR_API gErrorSize : public gErrorBaseHandler {
public:
  /**
   * 2D Size mismatch
   */
  gErrorSize(std::string classname,
             std::string funcname,
             std::string errname,
             T_Index rows1,
             T_Index cols1,
             T_Index rows2,
             T_Index cols2);
  /**
   * 1D Size mismatch
   */
  gErrorSize(std::string classname,
             std::string funcname,
             std::string errname,
             T_Index elem1,
             T_Index elem2);
  /**
   * 1D - 2D Size mismatch
   */
  gErrorSize(std::string classname,
             std::string funcname,
             std::string errname,
             T_Index rows1,
             T_Index cols1,
             T_Index elem2 );
}; // gErrorSize

/**
 * Slice resize error
 */
class GLINEAR_API gErrorSliceResize : public gErrorBaseHandler {
public:
  /**
   * Slice resize error
   */
  gErrorSliceResize(std::string classname,
                    std::string funcname);
}; // gErrorSliceResize

/**
 * Invalid argument error
 */
class GLINEAR_API gErrorInvalidArgument : public gErrorBaseHandler {
public:
  gErrorInvalidArgument(std::string classname,
                        std::string funcname,
                        std::string errstr);
}; // gErrorInvalidArgument

/**
 * Memory allocation failed
 */
class GLINEAR_API gErrorBadAlloc : public gErrorBaseHandler {
public:
  gErrorBadAlloc(std::string classname,
                 std::string funcname,
                 std::string errstr,
                 T_Index atom_size,
                 T_Index atoms);
}; // gErrorBadAlloc


GLINEAR_NAMESPACE_END

#endif // GERROR_H
