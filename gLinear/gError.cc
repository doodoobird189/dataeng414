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

char GERROR_CC_VERSION[] = "$Id$";

// gLinear headers
#include "gError.h"
#include "exceptions.h"

// standard headers
#include <string>  // string
#include <sstream>  // ostringstream
#include <iostream>  // endl, ends
#ifdef ABORT_NOT_THROW
#include <cstdlib>  // abort
#endif // ABORT_NOT_THROW

GLINEAR_NAMESPACE_BEGIN

// base handler
void gErrorBaseHandler::handle(const std::string& errstr,
             int) {
    std::ostringstream errMsg;
    errMsg << '\n' << "ERROR => " << errstr << std::endl;
    errMsg << __FILE__ << ":" << __LINE__ << std::endl;
#ifdef ABORT_NOT_THROW
    std::cerr << errMsg.str() << std::flush;
    abort();
#else
    throw GLinearError( errMsg.str() );
#endif // ABORT_NOT_THROW
}

// 2D index err
gErrorIndex::gErrorIndex(std::string classname,
       std::string funcname,
       T_Index maxrows,
       T_Index maxcols,
       T_Index row,
       T_Index col) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Index Error: Index out of bounds" << std::endl <<
     "  Allowed range: (0 .. " << maxrows << ",0 .. " << maxcols << ")" <<
     std::endl <<
     "  Actual values (" << row << "," << col << ")" << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// 1D index err
gErrorIndex::gErrorIndex(std::string classname,
       std::string funcname,
       std::string errname,
       T_Index maxval,
       T_Index val) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Index error: " << errname << std::endl <<
     "  Allowed range: 0 .. " << maxval << std::endl <<
     "  Actual values " << val << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// slice error
gErrorSlice::gErrorSlice(std::string classname,
       std::string funcname,
       std::string errname,
       T_Index maxindx,
       T_Index strt,
       T_Index end,
       T_Index stride) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Invalid stride: " << errname << std::endl <<
     "  Allowed range: 0 =< start =< end =< " << maxindx <<
     " stride > 0 " << std::endl <<
     "  Actual values: " << strt << " " << end << " " << stride << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// 1D and 2D subrange
gErrorSubrange::gErrorSubrange(std::string classname,
             std::string funcname,
             std::string errname,
             T_Index from,
             T_Index to,
             T_Index maxval) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Subrange Error: " << errname << std::endl <<
     "  Actual values: " << from << " " << to << " " << maxval << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// aliases exist
gErrorAliases::gErrorAliases(std::string classname,
           std::string funcname,
           std::string errname) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Aliases exist, operation not permitted: " << errname << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// shape error
gErrorShape::gErrorShape(std::string classname,
       std::string funcname,
       std::string errname,
       T_Index rows,
       T_Index cols) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Shape error: " << errname << std::endl <<
     "  [" <<
     rows << " x " << cols << "] unacceptable " << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// 2D size error
gErrorSize::gErrorSize(std::string classname,
           std::string funcname,
           std::string errname,
           T_Index rows1,
           T_Index cols1,
           T_Index rows2,
           T_Index cols2) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Size mismatch: " << errname << std::endl <<
     "  [" <<
     rows1 << " x " << cols1 << "] != [" <<
     rows2 << " x " << cols2 << "]" << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// 1D size error
gErrorSize::gErrorSize(std::string classname,
           std::string funcname,
           std::string errname,
           T_Index elem1,
           T_Index elem2) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Size mismatch: " << errname << std::endl <<
     "  [" << elem1 <<  "] != [" <<
     elem2 << "]" << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// 1D - 2D size error
gErrorSize::gErrorSize(std::string classname,
           std::string funcname,
           std::string errname,
           T_Index rows1,
           T_Index cols1,
           T_Index elem2 ) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Size mismatch: " << errname << std::endl <<
     "  [" <<
     rows1 << " x " << cols1 << "] != [" <<
     elem2 << "]" << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// Slice resize
gErrorSliceResize::gErrorSliceResize(std::string classname,
             std::string funcname) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Resize of slice type not allowed " << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// invalid argument
gErrorInvalidArgument::gErrorInvalidArgument(std::string classname,
               std::string funcname,
               std::string errstr) {
   std::ostringstream tmp;
   tmp <<
     "  In " << classname << "::" << funcname << std::endl <<
     "  Invalid argument: " << errstr << std::endl <<
     std::ends;
   handle( tmp.str() );
}

// memory allocation failed
gErrorBadAlloc::gErrorBadAlloc(std::string classname,
             std::string funcname,
             std::string errstr,
             T_Index atom_size,
             T_Index atoms) {
  std::ostringstream tmp;
  tmp <<
    "  In " << classname << "::" << funcname << std::endl <<
    "  Error: " << errstr << " allocating [ " << atoms << " of " << atom_size << " ] = " <<
    atoms * atom_size << " bytes" << std::endl <<
    std::ends;
  handle( tmp.str() );
}

GLINEAR_NAMESPACE_END
