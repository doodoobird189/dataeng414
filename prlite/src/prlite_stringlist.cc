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

char STRINGLIST_CC_VERSION[] = "$Id$";

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

// patrec headers
#include "prlite_stringlist.hpp"
#include "prlite_fnamepar.hpp"      // for changePathAndSuffix
#include "prlite_gendefs.hpp"        // for directory slash def
#include "prlite_error.hpp"

// standard headers
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>  // numeric_limits
#include <cctype>  // isspace
#include <new>  // placement new
#include <chrono>       // std::chrono::system_clock
#include <random>
//#include <iterator>


using namespace std;

//namespace prlite{

/************************** vector<string> streaming *************************/

// Insert STL vector of strings in stream
ostream& operator<<(ostream& os,
                    const vector<string>& vec) {
  unsigned N = vec.size();
  os << N << endl;
  for (unsigned i = 0; i < N; i++) {
    os << vec[i] << ' ';
  } // for i
  return os << endl;
} // operator<<

// Extract STL vector of strings from stream
istream& operator>>(istream& is,
                    vector<string>& vec) {
  size_t N;
  is >> N;
  vec.reserve(N);
  for (size_t i = 0; i < N; i++) {
    is >> vec[i];
  } // for i
  return is;
} // operator>>

/************************** StringList               *************************/

static const string TERMSTRING("\\\\\\");

//========================= Constructors             ========================

// Default constructor
StringList::StringList(int props)
    : slist(),
      listFile() {
  setProperties(props);
} // constructor

// Copy constructor
StringList::StringList(const StringList& p)
    : slist(p.slist),
      listFile(p.listFile) {
  // Make sure properties are consistent
  int props = p.getProperties();
  // (silently reset uniqueness flag)
  // This behaviour is useful when copying the StringList part of SearchableStringList,
  // which happens to always be 'unique' and would therefore raise a stink everytime
  if (props & 2) {
    props -= 2;
  } // if
  // The same applies to case insensitivity (can't be guaranteed from now on)
  if (props & 4) {
    props -= 4;
  } // if
  setProperties(props);
} // copy constructor

// Construct from file
StringList::StringList(const string& fileName,
                       int props)
    : slist(),
      listFile() {
  setProperties(props);
  loadFromFile(fileName);
} // constructor

// Construct from list of strings
StringList::StringList(const vector<string>& lst,
                       int props)
    : slist(),
      listFile() {
  setProperties(props);
  slist.reserve( lst.size() );
  // This ensures that each string is properly checked
  for (unsigned i = 0; i < lst.size(); i++) {
    add( lst[i] );
  } // for
} // constructor

// Destructor
StringList::~StringList() {
} // destructor

//========================= Operators and Conversion ========================

// Return immutable list element
const string& StringList::operator[](unsigned i) const {
  return slist[i];
} // operator[]

//========================= Exemplar Support         ========================

// Copy object shallowly
StringList* StringList::copy() const {
  return new StringList(*this);
} // copy

// Identify object
const string& StringList::isA() const {
  static const string CLASSNAME("StringList");
  return CLASSNAME;
} // isA

//========================= Inplace Configuration    ========================

// Configure via prompt
unsigned StringList::configure(unsigned cfgKey) {
  // Use default (unhelpful) helpstring
  return configureListOf("strings", cfgKey);
} // configure

// Configure via prompt with helpstring
unsigned StringList::configureListOf(const string& helpStr,
                                     unsigned cfgKey) {
  // Display help string
  cout << endl << "[01;35m" << "Configuring a " << isA() << "[00m";
  if (helpStr != "strings") {
    cout << " containing " << helpStr;
  } // if
  cout << endl;

  // Determine string capabilities
  cout << "Select the desired string properties:" << endl
       << "  (a) Do you want to allow the " << helpStr << " to contain whitespace (Y / N)? ";
  char ch;
  cin >> ch;
  bool allowSpaces = ( (ch == 'Y') || (ch == 'y') );
  int props = allowSpaces ? 1 : 0;
  // Searchable specific configures
  if (cfgKey == 1) {
    // Strings will be forced to be unique later on
    //    props += 2;
    // This question only makes sense with unique strings
    cout << "  (b) Are the " << helpStr << " case-sensitive (Y / N)? ";
    cin >> ch;
    props += (ch == 'N') || (ch == 'n') ? 4 : 0;
  } // if

  // Interactively determine list file
  bool fileNotFound = true;
  string fileName;
  do {
    cout << "Supply the name of the file containing the list of " << helpStr << " (. to create): ";
    cin >> fileName;

    // Create list from scratch
    if ( (fileName.size() == 1) && (fileName[0] == '.') ) {
      fileNotFound = false;
      vector<string> invec;
      string tmp;

      // Retrieve list of strings
      if (allowSpaces) {
        cout << "Enter the elements of the list one line at a time," << endl
             << "and terminate the process by entering a blank line:" << endl;
        cin.ignore(1, '\n');
        getString(cin, tmp);
        while (tmp.size() > 0) {
          invec.push_back(tmp);
          getString(cin, tmp);
        } // while
      } // if
      else {
        cout << "Enter the elements of the list separated by any form of whitespace," << endl
             << "and terminate the process by entering the string '" << TERMSTRING << "':" << endl;
        cin >> tmp;
        while (tmp != TERMSTRING) {
          invec.push_back(tmp);
          cin >> tmp;
        } // while
      } // else

      unsigned code = configure(invec, props);
      cout << "[01;34m" << "Finished configuring a " << isA() << "[00m" << endl;
      return code;
    } // if

    // Check file existence
    ifstream findfile(fileName.c_str(), ios::in);
    if (!findfile) {
      cout << "Sorry, cannot open file '" << fileName << "'" << endl;
    } // if
    else {
      findfile.close();
      fileNotFound = false;
    } // else
  } while (fileNotFound);

  unsigned code = configure(fileName, props);
  cout << endl << "[01;34m" << "Finished configuring a " << isA() << "[00m" << endl;
  return code;

} // configureListOf

// Configure from file
unsigned StringList::configure(const string& fileName,
                               int props,
                               unsigned) {
  // Destroy myself
  this->~StringList();
  // and do an inplace construction
  new(this) StringList(fileName, props);
  // Return success
  return 1;
} // configure

// Configure from list of strings
unsigned StringList::configure(const vector<string>& lst,
                               int props,
                               unsigned) {
  // Destroy myself
  this->~StringList();
  // and do an inplace construction
  new(this) StringList(lst, props);
  // Return success
  return 1;
} // configure

// Construct from text file
istream& StringList::txtRead(istream& file) {
  // Load string capabilities
  int props;
  file >> props;
  bool allowSpaces = ( (props % 2) == 1 );

  // Load list filename
  string fileName;
  file >> fileName;
  // If no list file specified
  if ( (fileName.size() == 1) && (fileName[0] == '-') ) {
    // Load list explicitly
    vector<string> lst;
    int numStrings;
    file >> numStrings;
    lst.resize(numStrings);
    if (allowSpaces) {
      file.ignore(1, '\n');
      for (int n = 0; n < numStrings; n++) {
        getString( file, lst[n] );
      } // for n
    } // if
    else {
      for (int n = 0; n < numStrings; n++) {
        file >> lst[n];
      } // for n
    } // else
    configure(lst, props);
  } // if
  else {
    configure(fileName, props);
  } // else
  return file;
} // txtRead

// Save to text file
ostream& StringList::txtWrite(ostream& file) const {
  // Save string capabilities
  file << endl << getProperties() << ' ';

  // Save list filename, if any
  if (listFile.size() > 0) {
    file << listFile << endl;
  } // if
  else {   // Else save list explicitly
    file << "- " << slist.size() << endl << *this << endl;
  } // else
  return file;
} // txtWrite

//========================= Required Members         ========================

// Number of elements in list
unsigned StringList::size() const {
  return slist.size();
} // size

// Check if empty
bool StringList::empty() const {
  return slist.empty();
} // empty

// Get properties
int StringList::getProperties() const {
  int props = 0;
  props += allowBlanks ? 1 : 0;
  props += unique ? 2 : 0;
  props += caseInsensitive ? 4 : 0;
  return props;
} // getProperties

// Clear list
void StringList::clear() {
  slist.clear();
  listFile = "";
} // clear

// Add string to list
void StringList::add(const string& str,
                     unsigned index) {
  // Make sure string is kosher
  checkProperties(str, index);
  // Append to list
  if ( index >= slist.size() ) {
    slist.push_back(str);
  } // if
  else {   // Else replace existing element
    slist[index] = str;
  } // else
} // add

// Add StringList to list
void StringList::add(const StringList& strlist,
                     unsigned index) {
  // Add strings one for one - necessary for checking
  for (unsigned n = 0; n < strlist.size(); n++) {
    add(strlist[n], index);
    // Advance index if inside list
    if ( index < numeric_limits<unsigned int>::max() ) {
      index++;
    } // if
  } // for n
} // add

// Remove string from list
void StringList::remove(unsigned index) {
  // Now shift strings up in list and then remove last guy
  for (unsigned i = index; i < slist.size() - 1; i++) {
    slist[i] = slist[i + 1];
  } // for i
  slist.pop_back();
} // remove

// Sort list
void StringList::sort() {
  std::sort( slist.begin(), slist.end() );
} // sort

void StringList::randomize() {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  //std::random_shuffle( slist.begin(), slist.end(), std::default_random_engine(seed));
   std::shuffle( slist.begin(), slist.end(), std::default_random_engine(seed));
} // randomize

// Load strings from specified file
void StringList::loadFromFile(const string& fileName) {
  // Open file
  ifstream file(fileName.c_str(), ios::in);
  PRLITE_CHECK(file, IOError, isA() << "::loadFromFile : cannot open file '" << fileName << "'");

  // Retrieve strings
  string tmp;
  if (allowBlanks) {
    getString(file, tmp);
    while (tmp.size() > 0) {
      add(tmp);
      getString(file, tmp);
    } // while
  } // if
  else {
    // Check for end of file
    while (file) {
      // Make sure tmp is empty
      tmp.erase();
      // Get next string
      file >> tmp;
      if (tmp.size() > 0) {
        add(tmp);
      } // if
    } // while
  } // else
  file.close();

  // Set list file name only on success
  listFile = fileName;
} // loadFromFile

// Save strings to specified file
void StringList::saveToFile(const string& fileName) const {
  // Open file
  ofstream file(fileName.c_str(), ios::out | ios::trunc);
  PRLITE_CHECK(file, IOError, isA() << "::saveToFile : cannot open file '" << fileName << "'");

  // Save list ALWAYS as newline-delimited, thereby compatible with both
  // allowBlanks and !allowBlanks
  for (unsigned n = 0; n < size(); n++) {
    file << slist[n] << endl;
  } // for n
  file.close();
} // saveToFile

// Print contents
void StringList::print(int depth) const {
  if (depth < 0) {
    return;
  } // if
  cout << isA();
  if ( listFile.size() ) {
    cout << " based on file '" << listFile << "'";
  } // if
  cout << " contains " << size() << " element(s)" << endl;

  if (depth < 1) {
    return;
  } // if
  cout << "String properties: ";
  if (allowBlanks) {
    cout << "may contain spaces, ";
  } // if
  else {
    cout << "no spaces, ";
  } // else
  if (unique) {
    cout << "unique, ";
  } // if
  else {
    cout << "non-unique, ";
  } // else
  if (caseInsensitive) {
    cout << "case-insensitive";
  } // if
  else {
    cout << "case-sensitive";
  } // else
  cout << endl;
  for (unsigned i = 0; i < slist.size(); i++) {
    cout << i << ": " << slist[i] << endl;
  } // for i
} // print

// Return data structure
const vector<string>& StringList::strings() const {
  return slist;
} // strings

// Return list filename
const string StringList::listFileName() const {
  return listFile;
} // listFileName

// Blank support
bool StringList::mayContainBlanks() const {
  return allowBlanks;
} // mayContainBlanks

// Uniqueness requirement
bool StringList::uniqueStrings() const {
  return unique;
} // uniqueStrings

// Case sensitivity
bool StringList::isCaseSensitive() const {
  return !caseInsensitive;
} // isCaseSensitive

//========================= Internal Members         ========================

// Set properties
void StringList::setProperties(int props) {
  allowBlanks = ( (props % 2) == 1 );
  // It's just too much effort to sustain a unique list without searchability
  unique = false;
  props >>= 1;
  if ( (props % 2) == 1 ) {
    cerr << "\nWARNING => " << isA() << "::setProperties :" << endl
         << "If you want a list with unique strings, use SearchableStringList instead" << endl;
  } // if
  props >>= 1;
  caseInsensitive = ( (props % 2) == 1 );
} // setProperties

// Check properties
void StringList::checkProperties(const string& str,
                                 unsigned) const {
  // Check if string illegally contains blanks
  // This is VERY important to ensure correct load / save of list
  if ( !allowBlanks && (str.find(' ') != string::npos) ) {
    PRLITE_THROW(StringListError,
          isA() << "::checkProperties for list '" << listFileName() << "':" << endl
          << "Cannot accept string '" << str << "' since it contains spaces");
  } // if
} // checkProperties

// Get string from stream, strictly delimited by a newline
void StringList::getString(istream& strm,
                           string& tmp) const {
  char ch;
  tmp.erase();

  // skip whitespace except '\n' up to first character
  do {
    if ( !strm.get(ch) ) {
      return;
    } // if
  } while ( ch != '\n' && isspace(ch) );
  // empty string so far
  if (ch == '\n') {
    return;
  } // if
  // otherwise first valid character hit
  tmp += ch;

  // read up to '\n', INCLUDING spaces (not tabs, etc. though)
  do {
    if ( !strm.get(ch) || (ch == '\n') ) {
      break;
    } // if
    if ( !isspace(ch) || (ch == ' ') ) {
      tmp += ch;
    } // if
  } while (ch != '\n');

  // chomp any final whitespace in string
  while ( (tmp.size() > 0) && isspace( *tmp.rbegin() ) ) {
    tmp.erase(tmp.size() - 1, 1);
  } // while
} // getString

//========================= Friends                  ========================

// Stream insertor
ostream& operator<<(ostream& os,
                    const StringList& lst) {
  // If strings contain blanks, delimit with '\n'
  if (lst.allowBlanks) {
    for (unsigned i = 0; i < lst.slist.size(); i++) {
      os << lst.slist[i] << endl;
    } // for i
  } // if
  else {  // Else delimit with ' '
    for (unsigned i = 0; i < lst.slist.size(); i++) {
      os << lst.slist[i] << ' ';
    } // for i
  } // else
  return os;
} // operator<<

// Stream extractor
istream& operator>>(istream& is,
                    StringList& lst) {
  vector<string> vec;
  is >> vec;
  lst.configure(vec, DEFAULT_SLIST_PROPS);
  return is;
} // operator>>

/*****************************************************************************/

// Modify path and suffix in filename list
void changePathAndSuffix(StringList& fileList,
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

  prlite::FileNameParse fparse;
  // Iterate through files in list
  for (unsigned i = 0; i < fileList.size(); i++) {
    // Construct new file name
    fparse.setFileName( fileList[i] );
    fileList.add( (keepPath ? fparse.path() : outPath) + fparse.main() + outSuffix, i );
  } // for i
} // changePathAndSuffix

//} // prlite
