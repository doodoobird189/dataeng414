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

char LOGGING_CC_VERSION[] = "$Id$";

/*
 * Authors    : Albert Visagie
 * Created on : 14 February 2005
 * Copyright  : University of Stellenbosch, all rights retained
 * Purpose    : Logging class.
 */

// patrec headers
#include "prlite_logging.hpp"
#include "prlite_error.hpp"
#include "prlite_stringops.hpp"  // split
#include "prlite_thread.hpp"

// standard headers
#include <ctime>  // ctime
#include <cstdlib>  // getenv

#ifdef __UNIX__
#include <unistd.h> // getpid
#else
// TODO getpid on other platforms?
#endif

using namespace std;

//namespace prlite{
/**
 * protect global things.
 */
static RecursiveMutex& __logger_mutex() {
  static RecursiveMutex __mutex;
  return __mutex;
}


/**
 * patrec's global logger.  by default it has no subscribers. applications
 * should subscribe to streams they want.
 */
Logger& getPatrecLogger() {
  // this is not thread safe: the first call should not be made in a
  // multithreaded environment.  Thereafter the logger protects
  // itself.
  static Logger __pLog;
  return __pLog;
}


LogWriter::~LogWriter() {

}


CerrLogWriter::CerrLogWriter(bool d) :
writeDate(d)
{

}

CerrLogWriter::~CerrLogWriter() {

}

void CerrLogWriter::logit(int evtType, const string& name, const LogEvent& evt) {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  string thetime;
  if (writeDate) {
    time_t t;
    time(&t);
    thetime = ctime(&t);
    thetime = thetime.substr(0,thetime.length()-1);
  }
  else {
    thetime = "";
  }
  cerr << name << "|";
  // note: if more than on stream is switched on, they will all appear
  // here.
  if (evtType & Logger::RUNNING_LOG) cerr << "LOGRUNNING";
  if (evtType & Logger::MORE_INFO) cerr << "LOGMOREINFO";
  if (evtType & Logger::DEBUG_LEVEL_1) cerr << "LOGDEBUG1";
  if (evtType & Logger::DEBUG_LEVEL_2) cerr << "LOGDEBUG2";
  if (evtType & Logger::ERROR_LOG) cerr << "LOGERROR";
  if (evtType & Logger::WARN) cerr << "LOGWARN";
  cerr << "|" << thetime << "|" << evt.message << endl;
  cerr << flush;
}


FileLogWriter::FileLogWriter(const string& filename, bool d, bool a) :
writeDate(d)
{
  if (a) {
    of.open(filename.c_str(), ios::binary | ios::app);
    PRLITE_CHECK(of, IOError,
          "Could not open " << filename << " in append mode for logging.");
  }
  else {
    of.open(filename.c_str(), ios::binary | ios::out);
    PRLITE_CHECK(of, IOError,
          "Could not create " << filename << " for logging.");
  }
  of.precision(12);
}

FileLogWriter::~FileLogWriter() {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  of.close();
}

void FileLogWriter::logit(int evtType, const string& name, const LogEvent& evt) {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  PRLITE_CHECK(of, IOError, "Cannot write to logfile!");
  string thetime;
  if (writeDate) {
    time_t t;
    time(&t);
    thetime = ctime(&t);
    thetime = thetime.substr(0,thetime.length()-1);
  }
  else {
    thetime = "";
  }
  of << name << "|";
  // note: if more than on stream is switched on, they will all appear
  // here.
  if (evtType & Logger::RUNNING_LOG) of << "LOGRUNNING";
  if (evtType & Logger::MORE_INFO) of << "LOGMOREINFO";
  if (evtType & Logger::DEBUG_LEVEL_1) of << "LOGDEBUG1";
  if (evtType & Logger::DEBUG_LEVEL_2) of << "LOGDEBUG2";
  if (evtType & Logger::ERROR_LOG) of << "LOGERROR";
  if (evtType & Logger::WARN) of << "LOGWARN";
  of << "|" << thetime << "|" << evt.message << endl;
  of << flush;
}


Logger::Logger() : mask(0) {

}

Logger::Logger(const Logger& l) :
subscribers(l.subscribers),
mask(l.mask)
{

}

Logger::~Logger() {

}

void Logger::subscribe(int evtClass, string streamName, RCPtr<LogWriter> sub) {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  SubscriberBlob blob;
  blob.name = streamName;
  blob.mask = evtClass;
  mask |= blob.mask;
  blob.sub = sub;
  subscribers.push_back(blob);
}

bool Logger::unsubscribe(int evtClass, string streamName) {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  bool ret = false;
  vector<SubscriberBlob>::iterator iter = subscribers.begin();
  while ( iter != subscribers.end() ) {
    if ( (iter->mask & evtClass) && (iter->name.find(streamName) == 0) ) {
      subscribers.erase(iter);
      ret = true;
      iter = subscribers.begin();
    }
    else {
      iter++;
    }
  }
  mask = 0;
  for ( ; iter!=subscribers.end(); iter++) {
    mask |= iter->mask;
  }
  return ret;
}

bool Logger::willLogit(int evtClass, const string& streamName) {
  RecursiveMutex::scoped_lock lock(__logger_mutex());
  bool ret = false;
  // quickquick test
  if (mask & evtClass) {
    // ok, test all
    vector<SubscriberBlob>::iterator iter = subscribers.begin();
    for ( ; iter!=subscribers.end(); iter++) {
      if ( (iter->mask & evtClass) &&
           ( (iter->name.find(streamName) != iter->name.npos) || (iter->name == "") ) ) {
        ret = true;
      }
    }
  } // if
  return ret;
}

void initLogging(const string& _procname, const string& _filename)
{
  volatile static bool initd = false;
  if (!initd) {
    const char* _value;
    int level = Logger::LEVEL1;
    string filename;
    string logLevelStr = "1";

    // get the logging level
    _value = getenv("LOGLEVEL");
    if (_value) {
      string logLevelStr = _value;
      if (logLevelStr == "0") {
        // don't subscribe anything
        return;
      }
      else if (logLevelStr == "1") {  level = Logger::LEVEL1;      }
      else if (logLevelStr == "2") {  level = Logger::LEVEL2;      }
      else if (logLevelStr == "3") {  level = Logger::LEVEL3;      }
      else if (logLevelStr == "4") {  level = Logger::LEVEL4;      }
      else                         {  level = Logger::LEVEL4;      }
    }

    // get the specific class that people wanted to subscribe to.
    _value = getenv("LOGNAMES");
    vector<string> logNames;
    if (_value) {
      string value = _value;
      prlite::split(value, ',', logNames);
    }

    // see if they wanted to log somewhere in particular.
    _value = getenv("LOGFILENAME");
    filename = _filename;
    if (_value) {
      filename = _value;
      if (filename == "") {
        filename = "-";
      }
    }

    RCPtr<LogWriter> writer;
    try {
      if (filename == "-") {
        writer = RCPtr<LogWriter>(new CerrLogWriter(true));    // write the date
      }
      else {
        writer = RCPtr<LogWriter>(new FileLogWriter(filename,
                                                    true,      // write the date
                                                    true));    // append to the existing file
      }
    }
    catch (exception& e) {
      cerr << "==========================================" << endl
           << "Exception when starting up logging system: " << endl
           << e.what() << endl;
    }

    PRLITE_ASSERT(writer, "Something bad is wrong with the logging initialization code.");

    if (logNames.size() == 0) {
      // nothing given, subscribe to everything.
      PRLITE_pLog.subscribe(level, "", writer);
    }
    else {
      // some given in the environment variable, subscribe to them.
      vector<string>::iterator iname = logNames.begin();
      for ( ; iname!=logNames.end(); iname++) {
        PRLITE_pLog.subscribe(level, *iname, writer);
      }
    }

    // let them know we're up
#ifdef DEBUG_MODE
    static const string name = "PatrecII_DEBUG";
#else
    static const string name = "PatrecII";
#endif

    ostringstream oss;
#ifdef __UNIX__
    oss << "Running " << _procname << " (" << getpid() << "). Logging system started. "
#else
    oss << "Running " << _procname << ". Logging system started. "
#endif
        << "LOGFILENAME=\"" << filename << "\", "
        << "LOGLEVEL=" << logLevelStr << ", "
        << "LOGNAMES=\"" << prlite::join(",", logNames.begin(), logNames.end()) << "\"";
    writer->logit(Logger::RUNNING_LOG, name,
                  LogEvent(oss.str()));

    initd = true;
  } // if (!initd)
}

//} // namespace prlite
