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

char MY_ANY_CC_VERSION[] = "$Id$";

/*******************************************************************************

   AUTHOR:  Albert Visagie (Copyright University of Stellenbosch, all rights reserved.)
   DATE:    11 August 2005
   PURPOSE: utlities for Any class.

*******************************************************************************/

// patrec includes
#include "prlite_my_any.hpp"

// standard includes
#include <string>

using namespace std;

//namespace prlite{

/**
 * useful little helper: constructs an Any according to the type string.
 * This should probably become a standard utility for use with Any objects.
 */
Any mkAny(const string& value, const string& type) {
  istringstream iss(value);
  Any ret;
  if (type == "int") {
    int tmp;
    iss >> tmp;
    ret = tmp;
  }
  else if (type == "float") {
    float tmp;
    iss >> tmp;
    ret = tmp;
  }
  else if (type == "string") {
    ret = value;
  }
  else {
    PRLITE_FAIL("Only \"int\", \"float\" and \"string\" recognised as types!" << endl
         << "Got \"" << type << "\".");
  }
  PRLITE_ASSERT(iss, "Could not convert \"" << value << "\" to " << type << "!");
  return ret;
}

string anyType(const Any& any)
{
  if (any.type() == typeid(int)) {
    return "int";
  }
  else if (any.type() == typeid(float)) {
    return "float";
  }
  else if (any.type() == typeid(string)) {
    return "string";
  }
  else {
    PRLITE_FAIL("Only \"int\", \"float\" and \"string\" reconised as types!");
  }
  // To keep the compiler happy
  return "";
}

//} // namespace prlite
