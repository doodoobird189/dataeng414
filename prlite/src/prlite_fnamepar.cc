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

char FNAMEPAR_CC_VERSION[] = "$Id$";

 /*****************************************************************************
 *               |                                                            *
 * Author        | Johan du Preez, Ludwig Schwardt                            *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 1995/09/14                                                 *
 * Description   | A facility to split filenames in various parts             *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

// patrec headers
#include "prlite_fnamepar.hpp"
#include "prlite_gendefs.hpp"        // for directory slash def

using namespace std;

/************************** FileNameParse            *************************/

namespace prlite{

FileNameParse::FileNameParse() {
} // default constructor

FileNameParse::FileNameParse(const string& fname) {
  setFileName(fname);
} // constructor

// Set file name to parse
void FileNameParse::setFileName(const string& fname) {
  string fullName(fname);
  int punt = fullName.rfind('.', fullName.length() - 1);
  int slash;
  if (punt < 0) {
    punt = fullName.length();
    slash = fullName.rfind(DIRSLASH, fullName.length() - 1);
  } // if
  else {
    slash = fullName.rfind(DIRSLASH, punt);
  } // else
  pth = fullName.substr(0, slash + 1);
  int ssLen = punt - 1 - slash;
  mn = fullName.substr(slash + 1, ssLen);
  ssLen = fullName.length() - punt + 1;
  sfx = fullName.substr(punt, ssLen);
} // setFileName

string FileNameParse::path() const {
  return pth;
} // path

string FileNameParse::main() const {
  return mn;
} // main

string FileNameParse::suffix() const {
  return sfx;
} // suffix

/*****************************************************************************/

// Modify path and suffix of filename
extern void changePathAndSuffix(string& fileName,
        const string& newSuffix,
                                const string& newPath) {
  // Check if path should remain the same
  bool keepPath = ( (newPath.size() == 1) && (newPath[0] == '-') );

  // Append slash to path if it's not there (more userfriendly that way)
  string outPath = newPath;
  if ( outPath.size() && (outPath[outPath.size() - 1] != DIRSLASH) ) {
    outPath += DIRSLASH;
  } // if
  // Prepend . to extension if it's not there (more userfriendly that way)
  string outSuffix = newSuffix;
  if ( outSuffix.size() && (outSuffix[0] != '.') ) {
    outSuffix = '.' + outSuffix;
  } // if

  // Construct new file name
  FileNameParse fparse;
  fparse.setFileName(fileName);
  fileName = (keepPath ? fparse.path() : outPath) + fparse.main() + outSuffix;
} // changePathAndSuffix


static string _join(const string& thing1, const string& thing2)
{
    if ((thing1.rfind('\\') != thing1.length()-1) &&
  (thing1.rfind('/') != thing1.length()-1)) {
        return thing1 + DIRSLASH + thing2;
    }
    else {
        // there is a trailing slash on thing1
        return thing1 + thing2;
    }
}


string pathjoin( const string& thing1,
     const string& thing2,
     const string& thing3,
     const string& thing4,
     const string& thing5)
{
    string ret = thing1;
    ret = _join(ret, thing2);
    if (thing3 != "")
         ret = _join(ret, thing3);
    if (thing4 != "")
         ret = _join(ret, thing4);
    if (thing5 != "")
         ret = _join(ret, thing5);
    return ret;
}

} // namespace prlite
