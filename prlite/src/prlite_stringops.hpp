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

#ifndef PRLITE_STRINGOPS_HPP
#define PRLITE_STRINGOPS_HPP

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14/Nov/99
         PURPOSE: Operations on strings

*******************************************************************************/

// patrec headers
#include "prlite_gendefs.hpp"  // PRLITE_API

// standard headers
#include <iosfwd>  // istream, ostream
#include <string>  // string
#include <vector>  // vector
#include <sstream>  // istringstream
#include <limits>  // numeric_limits

namespace prlite{

PRLITE_API
extern std::ostream& write( std::ostream& file,
                            const std::string& st );
PRLITE_API
extern std::istream& read( std::istream& file,
                           std::string& st,
                           int maxlen = std::numeric_limits<int>::max());
PRLITE_API extern std::string lower( const std::string& str );
PRLITE_API extern std::string upper( const std::string& str );

inline std::string operator*(const std::string& str, int ii) {
  std::string ret;
  ret.reserve((str.length()+1) * ii);
  while (ii--) { ret += str; }
  return ret;
}

/**
 * strip whitespace from the start and end of a string.
 */
PRLITE_API std::string strip(const std::string& str);

/**
 * Read all the lines in file into a vector.
 */
PRLITE_API std::vector<std::string> readlines(std::istream& is);

PRLITE_API
unsigned int numcharinstring(const std::string& str,
                             const char c);


/**
 * Split a string on a particular delimiter character.
 */
PRLITE_API
void split(const std::string& s,
           const char c,
           std::vector<std::string>& v);

// splits on white-space
PRLITE_API
void split(const std::string& s,
           std::vector<std::string>& v);

/**
 * Split on a string rather than a character. for example usage, see the
 * testcode in prlite_stringops.cc
 */
template<typename iterator>
     void split(const std::string& s,
    const std::string& delimiter,
    iterator iter,
    bool doStrip = false) {
       std::string::size_type pos( s.find(delimiter) );
       if (pos == std::string::npos) {
   *iter = doStrip ? strip(s) : s;
   ++iter;
       } // if
       else {
   std::string::size_type len( delimiter.length() );
   *iter = doStrip ? strip(s.substr(0, pos)) : s.substr(0, pos);
   ++iter;
   split(s.substr(pos + len, s.size() - pos - len), delimiter, iter, doStrip);
       } // else
     } // split

/**
 * split on whitespace
 */
template<typename iterator>
     void split(const std::string& s,
    iterator iter) {
       std::istringstream iss( s );
       std::string tmp;
       iss >> tmp;
       while (iss) {
   *iter = tmp;
   ++iter;
   iss >> tmp;
       }
     } // split


/**
 * helper: join some strings into one string using filler string.
 */
template<typename Iter>
std::string join(const std::string& filler, Iter ii, const Iter& end) {
  if (ii==end) {
    return "";
  }
  std::string ret = *ii;
  for (++ii; ii!= end; ii++) {
    ret += filler + *ii;
  }
  return ret;
}

} // namespace prlite

#endif // PRLITE_STRINGOPS_HPP
