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

#ifndef PRLITE_FNAMEPAR_HPP
#define PRLITE_FNAMEPAR_HPP

// patrec headers
#include "prlite_gendefs.hpp"  // PRLITE_API

// standard headers
#include <string>  // string


namespace prlite{

/**************************** FileNameParse         **************************/

/**
 * A facility to split filenames in various parts, namely:
 * a path (everything and including the final '/',
 * a main part (everything upto and excluding the final '.',
 * and a suffix (everything following, and including the final '.')
 */

class PRLITE_API FileNameParse {

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
public:

  /**
   * Default constructor.
   */
  FileNameParse();

  /**
   * Initialize with file name to parse.
   */
  FileNameParse(const std::string& fname);

  //========================= Operators and Conversions =========================

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

  //=========================== Inplace Configuration ===========================

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================

  //======================= Non-virtual Member Functions ========================
public:

  /**
   * Set file name to parse.
   */
  void setFileName(const std::string& fname);

  /**
   * Return file path - everything upto and including final slash.
   */
  std::string path() const;

  /**
   * Return main part of file name - the bit between the final slash and the
   * final '.'
   */
  std::string main() const;

  /**
   * Return file extension - everything after and including final '.').
   */
  std::string suffix() const;

  //================================== Friends ==================================

  //=============================== Data Members ================================
private:

  /**
   * Path/name/extension strings.
   */
  std::string pth, mn, sfx;

  //============================ DEPRECATED Members =============================


}; // FileNameParse


/*****************************************************************************/

/**
 * Modify path and suffix of filename.
 * @param fileName  Filename to be changed
 * @param newSuffix New suffix (mandatory)
 * @param newPath   New path ("-" if to be kept the same)
 */
PRLITE_API
extern void changePathAndSuffix( std::string& fileName,
                                 const std::string& newSuffix,
                                 const std::string& newPath = "-" );

/**
 * Concatenate paths. The last one does not get a trailing DIRSLASH. For
 * example, If you give it foo and bar you get foo/bar, and if you give it
 * foo/ and bar, you get foo/bar again.  foo and bar/ will yield foo/bar/,
 * so take care.
 *
 * Surrounding space is not stripped.  Keep your nose clean.
 */
PRLITE_API
extern std::string pathjoin( const std::string& thing1,
                             const std::string& thing2,
                             const std::string& thing3 = "",
                             const std::string& thing4 = "",
                             const std::string& thing5 = "");

} // namespace prlite

#endif // PRLITE_FNAMEPAR_HPP
