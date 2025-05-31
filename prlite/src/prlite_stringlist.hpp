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
 * Author        | Ludwig Schwardt                                            *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 2001/08/26                                                 *
 * Description   | A simple vector<string> class loadable from a file         *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

#ifndef STRINGLIST_HPP
#define STRINGLIST_HPP

// patrec headers
#include "prlite_error.hpp"  // DEFINE_EXCEPTION
#include "prlite_gendefs.hpp"  // PRLITE_API

// standard headers
#include <string>  // string
#include <vector>  // vector
#include <iosfwd>  // istream, ostream
#include <limits>  // numeric_limits
#include <ciso646>  // and etc

 //namespace prlite{

/// Default string properties == no blanks, non-unique, case sensitive
static const int DEFAULT_SLIST_PROPS(0);

PRLITE_DEFINE_EXCEPTION(StringListError, PatRecError);

/**************************** StringList            **************************/

/**
 * This class provides a list of strings that can be loaded from a
 * file.  In addition, the list of strings can also be directly
 * entered from stdin.  Unlike FileList, this class doesn't assume the
 * strings to have a special structure, and they can therefore be used
 * for filenames, ID strings, tags, whatever. In addition, the class
 * has many standard members, allowing easy integration within other
 * objects. The strings are stored in an STL vector for easy access,
 * as well as easy insertion/deletion compared to gLinear Vectors.
 * The price to be paid (obligatory deepness) doesn't seem as
 * compelling.  The strings are also allowed to have certain
 * properties, such as uniqueness, case (in)sensitivity and the
 * inclusion of spaces. The first two options are forced false since
 * it overlaps too much with the concept of a
 * SearchableStringList. The last option has ramifications for the
 * loading and saving of strings.
 *
 * TODO: - new abstract base class that only handles properties, get and virtuals
 *       - new Searchable that doesn't contain vector<string>, only map
 *       - new BidirList that contains vector + map - multiple inheritance?
 *       - or not (too convoluted)?
 *
 * @author Ludwig Schwardt
 * <schwardt@dsp.sun.ac.za>
 * Copyright: University of Stellenbosch, all rights retained.
 */

class PRLITE_API StringList {

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
public:

  /**
   * Default constructor.
   * @param props  Flags indicating string properties
   */
  StringList(int props = DEFAULT_SLIST_PROPS);

  /**
   * Copy constructor.
   */
  StringList(const StringList& p);

  /**
   * Construct from file.
   * Make sure its explicit, otherwise StringLists cannot be used where
   * strings are expected (such as in the add() member).
   * @param fileName  Name of file containing list
   * @param props  Flags indicating string properties, decoded as:
   *               +1 - Strings may contain spaces
   *               +2 - Strings must be unique
   *               +4 - Strings are case insensitive
   */
  explicit StringList(const std::string& fileName,
                      int props = DEFAULT_SLIST_PROPS);

  /**
   * Construct from vector of strings.
   * @param lst  Vector of strings
   * @param props  Flags indicating string properties
   */
  StringList(const std::vector<std::string>& lst,
             int props = DEFAULT_SLIST_PROPS);

  /**
   * Destructor.
   */
  virtual ~StringList();

  //========================= Operators and Conversions =========================
public:

  /**
   * Return specific list entry (immutable).
   * No bitching if out of bounds (compile with index checking enabled instead
   * if this is desired). It is not desirable to return a non-const ref, since
   * this screws up later down the line when string lists become searchable.
   * Modification of the list is controlled via explicit members to ensure the
   * string properties are not violated.
   * @return  The i'th string (immutable) from the list.
   */
  const std::string& operator[](unsigned i) const;

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

  //=========================== Inplace Configuration ===========================
public:

  /**
   * Inplace configuration from prompt.
   * @param cfgKey  Configuration option (0 for StringList, 1 for Searchable-)
   */
  virtual unsigned configure(unsigned cfgKey = 0);

  /**
   * Inplace configuration from prompt, with help string.
   * This is a bit more user-friendly, as it reminds the user what the list contains,
   * instead of just stating "strings". Without it, configuration questions can appear
   * very obscure. A special member name is required, to avoid ambiguity with
   * configure(fileName).
   * @param helpStr This string reminds the user what the list contains during questions
   * @param cfgKey  Configuration option (0 for StringList, 1 for Searchable-)
   */
  virtual unsigned configureListOf(const std::string& helpStr,
                                   unsigned cfgKey = 0);

  /**
   * Inplace configuration from file.
   * @param fileName  Name of file containing list
   * @param props  Flags indicating string properties
   */
  virtual unsigned configure(const std::string& fileName,
                             int props = DEFAULT_SLIST_PROPS,
                             unsigned cfgKey = 0);

  /**
   * Inplace configuration from list of strings.
   * @param lst  Vector of strings
   * @param props  Flags indicating string properties
   */
  virtual unsigned configure(const std::vector<std::string>& lst,
                             int props = DEFAULT_SLIST_PROPS,
                             unsigned cfgKey = 0);

  /**
   * Construct object from text file.
   */
  virtual std::istream& txtRead(std::istream& file);

  /**
   * Save object to text file.
   */
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================
public:


  /**
   * Makes shallow copy of object.
   */
  virtual StringList* copy() const;

  /**
   * Identify object.
   */
  virtual const std::string& isA() const;

  /**
   * Discard all strings.
   */
  virtual void clear();

  /**
   * Add a string to list.
   * This can grow the list by appending strings past the end, as well as
   * replace strings at specific indices within the list. It is basically a
   * add/replace combination, pooled together for simplicity's sake.
   * @param str  String to add (default = empty)
   * @param index  Position in list where string will be placed (default = append to list)
   */
  virtual void add( const std::string& str = "",
                    unsigned index = std::numeric_limits<unsigned int>::max() );

  /**
   * Remove string from list.
   * This shrinks the list by eliminating the desired position. All strings at
   * later indices move one place up, and the list becomes one shorter. By default
   * it removes strings from the end of the list.
   * @param index  Position of string to be removed (default = last in list)
   */
  virtual void remove( unsigned index = std::numeric_limits<unsigned int>::max() );

  /**
   * Sort strings alphabetically (lexicographically).
   * This does an inplace sorting of the strings (alphabetically, capitals first).
   */
  virtual void sort();

  /**
   * Randomize strings inplace
   */
  virtual void randomize();

  /**
   * Load strings from specified file.
   * Similar to configure member, kept for symmetry with saveToFile.
   * @param  List file name
   */
  virtual void loadFromFile(const std::string& fileName);

  /**
   * Print information to specified depth.
   * @param depth  Depth of detail required.
   */
  virtual void print(int depth = 0) const;

protected:

  /**
   * Decode and set string properties.
   * Set string flags based on encoded integer. Just add the appropriate
   * values below for all set flags to create the property integer (if
   * a flag is not set, ignore the value). This kept protected so that
   * properties can only be set during construction, and NOT during the
   * life of a StringList object.
   * @param props  Flags indicating string properties, decoded as:
   *               +1 - Strings may contain spaces
   *               +2 - Strings must be unique
   *               +4 - Strings are case insensitive
   */
  virtual void setProperties(int props);

  /**
   * Check string properties.
   * This will currently abort if you try to enter invalid strings into a list.
   * String properties are therefore seen as preconditions and not to be taken
   * lightly.
   * @param str  String to be checked
   * @param pos  Position where string is to be added
   */
  virtual void checkProperties( const std::string& str,
                                unsigned pos = std::numeric_limits<unsigned int>::max() ) const;

  //======================= Non-virtual Member Functions ========================
public:

  /**
   * Number of strings in the list.
   * Return Number of strings in the list.
   */
  unsigned size() const;

  /**
   * Check if list is empty.
   * Return True if list is empty.
   */
  bool empty() const;

  /**
   * Get encoded string properties.
   * This returns the string property flags encoded within the least significant
   * bits of an integer.
   */
  int getProperties() const;

  /**
   * Add a string list to the list.
   * This applies add(string, index) to the elements of a specified string list.
   * Each string is therefore fully checked before insertion. If the index points
   * within the list, the second list overwrites the portion of the list starting
   * at the specified index. Otherwise, the second list is appended to the list.
   * @param strlist  StringList to add
   * @param index  Position in list where strlist will be placed (default = append to list)
   */
  void add( const StringList& strlist,
            unsigned index = std::numeric_limits<unsigned int>::max() );

  /**
   * Save strings to specified file.
   * Allows backup of string list.
   * @param  List file name
   */
  void saveToFile(const std::string& fileName) const;

  /**
   * Return immutable version of internal string structure.
   */
  const std::vector<std::string>& strings() const;

  /**
   * Return list filename, if any applicable.
   */
  const std::string listFileName() const;

  /**
   * Indicate support for strings containing blanks.
   * @return  True if blanks supported
   */
  bool mayContainBlanks() const;

  /**
   * Indicate whether strings have to be unique.
   * @return  True if strings are required to be unique
   */
  bool uniqueStrings() const;

  /**
   * Indicate whether strings are case sensitive.
   * @return  True if case sensitive
   */
  bool isCaseSensitive() const;

protected:

  /**
   * Get string from stream.
   */
  void getString(std::istream& strm,
                 std::string& tmp) const;

  //================================== Friends ==================================
public:

  /**
   * Stream insertor.
   * This is used to provide a simple display of the contents of the string list.
   * If the strings contain no blanks, they are printed with spaces as delimiters;
   * otherwise newline delimiters are used. No mention is made of the listfile.
   */
  friend std::ostream& operator<<(std::ostream& os,
                                  const StringList& lst);

  /**
   * Stream extractor.
   * This provides very rudimentary input to the string list (simply here to
   * mirror the insertor above). No blanks are allowed, as this merely thunks
   * directly through to the vector >> underneath. I don't think there is a
   * good reason to use this...
   */
  friend std::istream& operator>>(std::istream& is,
                                  StringList& lst);

  //=============================== Data Members ================================
protected:

  /**
   * String list.
   */
  std::vector<std::string> slist;

  /**
   * List filename, if any applicable.
   */
  std::string listFile;

  /**
   * Flag indicating whether strings are allowed to contain spaces.
   * If this is true, strings in the list are allowed to contain spaces,
   * and are strictly delimited by newlines. The alternative is strings that
   * must consist of non-whitespace only, and are delimited by any whitespace
   * characters.
   */
  bool allowBlanks;

  /**
   * Flag indicating whether strings have to be unique.
   * If this is true, the list is not allowed to contain duplicates.
   * This will always be false, since a
   * StringList + uniqueness = SearchableStringList.
   */
  bool unique;

  /**
   * Flag indicating whether strings are case sensitive.
   * If this is true and strings are also unique, the list may not contain
   * two strings that differ only in case. This will always be false, since
   * StringLists are non-unique.
   */
  bool caseInsensitive;

  //============================ DEPRECATED Members =============================


}; // StringList


/*****************************************************************************/

/**
 * Send STL vector of strings to stream.
 */
std::ostream& operator<<(std::ostream& os,
                         const std::vector<std::string>& vec);

/**
 * Get STL vector of strings from stream.
 */
std::istream& operator>>(std::istream& is,
                         std::vector<std::string>& vec);

/**
 * Change the path and suffix in a filename list.
 * This modifies a filename list by replacing the suffix and optionally the
 * path as well. As a further attempt at user-friendliness it appends a slash
 * to the new path if it does not contain one at the end. It assumes the
 * StringList contains file names.
 * @param fileList  List of file names
 * @param newSuffix  New suffix
 * @param newPath  New path (default = '-' indicating the path is unchanged)
 */
extern void changePathAndSuffix( StringList& fileList,
                                 const std::string& newSuffix,
                                 const std::string& newPath = "-" );


//} // prlite

#endif // STRINGLIST_HPP
