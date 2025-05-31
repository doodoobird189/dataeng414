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

char STRINGOPS_CC_VERSION[] = "$Id$";

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14 Nov 99
         PURPOSE: Operations on strings

*******************************************************************************/

// patrec headers
#include "prlite_stringops.hpp"  // lower, upper, write, read
#include "prlite_error.hpp"
#include "prlite_endian.hpp"

// standard headers
#include <string>  // string
#include <algorithm>  // transform
#include <cctype>  // tolower, toupper, isspace
#include <iostream>  // ostream, istream

using namespace std;

namespace prlite{

string lower(const string& str) {
  string tmp(str);
  for (unsigned i = 0; i < tmp.length(); i++) {
    tmp[i] = static_cast<char>(tolower( tmp[i] ));
  } // for
  // transform(str.begin(), str.end(), tmp.begin(), tolower);
  return tmp;
} // lower

string upper(const string& str) {
  string tmp(str);
  for (unsigned i = 0; i < tmp.length(); i++) {
    tmp[i] = static_cast<char>(toupper( tmp[i] ));
  } // for
  // transform(str.begin(), str.end(), tmp.begin(), toupper);
  return tmp;
} // upper

ostream& write( ostream& os,
                const string& st ) {
  int len( toLittleEndian(st.length()) );
  os.write( reinterpret_cast<char *>(&len), sizeof(len) );
  os.write( reinterpret_cast<const char *>( st.c_str() ), len * sizeof(char) );
  return os;
} // write

istream& read( istream& is,
               string& st,
         int maxlen) {
  int len;
  is.read( reinterpret_cast<char *>(&len), sizeof(len) );
  if (is.bad()) { return is; }
  // check the number
  len = PRLITE_fromLittleEndian(len);
  PRLITE_CHECK(len >= 0, IOError,
        "Bad length " << len << " when reading raw string from file.");
  PRLITE_CHECK(len <= maxlen, IOError,
        "Bad length " << len << " when reading raw string from file.");

  char* ch = new char[len + 1];
  is.read( reinterpret_cast<char *>(ch), len * sizeof(char) );
  if (is.bad()) { delete[] ch; return is; }
  ch[len] = '\0';
  st = ch;

  delete[] ch;
  return is;
} // read

unsigned int numcharinstring(const string& str,
                             const char c) {
  unsigned int count(0);
  for (string::const_iterator i( str.begin() ); i != str.end(); i++) {
    if (*i == c) {
      count++;
    } // if
  } // for
  return count;
} // numcharinstring

void split(const string& s,
           const char c,
           vector<string>& v) {
  string::size_type p( s.find(c) );
  if (p == string::npos) {
    v.push_back(s);
  } // if
  else {
    v.push_back( s.substr(0, p) );
    split(s.substr(p + 1, s.size() - p - 1), c, v);
  } // else
} // split

void split(const string& s,
           vector<string>& v) {
  istringstream iss( s );
  string tmp;
  iss >> tmp;
  while (iss) {
    v.push_back(tmp);
    iss >> tmp;
  } // while
} // split


string strip(const string& s) {
  string::size_type p( s.find_last_not_of("\t \n\r") );
  if (p == string::npos) {
    return "";
  } // if
  else {
    string s2( s.substr(0, p + 1) );
    p = s.find_first_not_of("\t ");
    if (p == string::npos) {
      return "";
    } // if
    else {
      return s2.substr(p);
    } // else
  } // else
} // strip

/*
string strip(const string& str) {
  cout << __FILE__ << __LINE__ << endl;
  const unsigned int ssize( str.size() );

  // seek the beginning from the beginning
  unsigned int begin(0);
  while ( (begin < ssize) && isspace( str[begin] ) ) {
    begin++;
  } // while

  if (begin == ssize) {
    return "";
  } // if

  // seek the end from the end
  int end = (ssize - 1);
  while (end >= 0 && str[end] == '\n' || str[end] == '\r' || str[end] == (unsigned char)(13) ) {
    end--;
  } // if
  while ( (end >= 0) && isspace( str[end] ) ) {
    end--;
  } // while

  if (end < 0) {
    return "";
  } // if

  return str.substr(begin, end - begin + 1);
}
*/

/**
 * Read all the lines in file into a vector.
 */
vector<string> readlines(istream& is) {
  string tmp;
  vector<string> ret;

  if (is) {
    getline(is, tmp);
  }

  while (is) {
    ret.push_back(tmp);
    getline(is, tmp);
  }

  return ret;
} // readlines

} // namespace prlite
