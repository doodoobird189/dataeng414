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

/*
 * Authors    : Albert Visagie
 * Created on : 14 February 2005
 * Copyright  : University of Stellenbosch, all rights retained
 * Purpose    : Logging class.
 */

#ifndef PRLITE_LOGGING_HPP
#define PRLITE_LOGGING_HPP

// patrec headers
#include "prlite_rcptr.hpp"
#include "prlite_gendefs.hpp"

// standard headers
#include <string>
#include <vector>
#include <fstream>

 //namespace prlite{

/**
 * Log Event class.  the tags field could contain things like vectors, or
 * debugging info or whatever.  The LogWriters are responsible to dealing
 * with this.
 */
struct LogEvent {
  /// a normal log message.
  std::string message;

  /// in case logs should ever contain more, uncomment this.
  /// std::map<std::string, Any> tags;

  LogEvent() { }
  LogEvent(const std::string& msg) : message(msg) { }
};

/**
 * Subscribers to log streams must support this interface.
 */
class PRLITE_API LogWriter {
public:
  virtual ~LogWriter();
  virtual void logit(int evtType, const std::string& name, const LogEvent& evt) = 0;
};

/**
 * A simple one that writes LogEvent::name to cerr.
 */
class PRLITE_API CerrLogWriter : public LogWriter {
  bool writeDate;
public:
  CerrLogWriter(bool writeDate = true);
  virtual ~CerrLogWriter();
  virtual void logit(int evtType, const std::string& name, const LogEvent& evt);
};

/**
 * Another simple one that writes LogEvent::name to a file.
 */
class PRLITE_API FileLogWriter : public LogWriter {
  bool writeDate;
  std::ofstream of;
public:
  FileLogWriter(const std::string& filename, bool writeDate = true, bool append = true);
  virtual ~FileLogWriter();
  virtual void logit(int evtType, const std::string& name, const LogEvent& evt);
};

/**
 * Patrec logging interface.
 *
 * The proposed log event types are:
 *
 * RUNNING_LOG is simple notifiers like "server started", "server
 * shutdown".
 *
 * MORE_INFO could be things like "Synthesise Utterance Request: 'my
 * telephone number is 123231'"
 *
 * ERROR_LOG for all errors, would probably be used most often by exception
 * classes?
 *
 * DEBUG level 1 is for messages that are useful to getting the system
 * up and running.  things like file not found etc.  warnings etc.
 *
 * DEBUG level 2 is for messages that are pretty much only useful to a
 * programmer with the source code.  this could contain things like
 * linenumbers and source file name etc.
 *
 * For the moment patrec all patrec objects should call
 *
 * pLog.moreInfo(isA(), "message");
 *
 * or
 *
 * pLog.debug1(isA(), "message");
 *
 * or
 *
 * pLog.debug2(isA(), "message");
 *
 * or
 *
 * pLog.warn(isA(), "message");
 *
 * Where exceptions are caught, but not handled the following should be used.
 *
 * pLog.error(isA(), "message")
 *
 * something should always subscribe to WARN and ERROR_LOG
 * pLog.subscribe(Logger::WARN | Logger::ERROR_LOG, "", writer);
 *
 * at some stage CAUTION will be changed to reflect this.
 */
class PRLITE_API Logger {

public:
  /**
   * This defines the bitmap for the five different event types, These
   * should be used in calls to logit or it's shortcuts.
   */
  typedef enum { RUNNING_LOG = 1<<0,
   MORE_INFO = 1<<1,
   WARN = 1<<2,
   DEBUG_LEVEL_1 = 1<<3,
   DEBUG_LEVEL_2 = 1<<4,
   ERROR_LOG = 1<<5
  } Logging_Streams;

  /**
   * This predefines certain levels in terms of the various event types
   * that are enabled in those levels. These should only be used in calls
   * to Logger::subscribe.
   */
  typedef enum { LEVEL0 = 0, // no logging
   LEVEL1 = RUNNING_LOG | WARN | ERROR_LOG,
   LEVEL2 = RUNNING_LOG | WARN | MORE_INFO | ERROR_LOG,
   LEVEL3 = RUNNING_LOG | WARN | MORE_INFO | DEBUG_LEVEL_1 | ERROR_LOG,
   LEVEL4 = RUNNING_LOG | WARN | MORE_INFO | DEBUG_LEVEL_1 | DEBUG_LEVEL_2 | ERROR_LOG,
   // levels, but without ERROR_LOG
   LEVEL1NE = RUNNING_LOG | WARN ,
   LEVEL2NE = RUNNING_LOG | WARN  | MORE_INFO,
   LEVEL3NE = RUNNING_LOG | WARN  | MORE_INFO | DEBUG_LEVEL_1,
   LEVEL4NE = RUNNING_LOG | WARN  | MORE_INFO | DEBUG_LEVEL_1 | DEBUG_LEVEL_2
  } Logging_Groups;

public:

  /**
   * default constructor.  this results in a logger with no subscribers,
   * and hence no output will be sent anywhere.
   */
  Logger();

  /// copy
  Logger(const Logger& l);

  /// die
  ~Logger();

  /**
   * subscribe a log writer to a set of streams.  The same LogWriter can be
   * subscribed multiple times.  making streamName = "" subscibes the
   * writer to all events with appropriate evtTypes.
   */
  void subscribe(int evtClass, std::string streamName, RCPtr<LogWriter> sub);

  /**
   * unsubscribe a LogWriter object.  returns true if one or more
   * LogWriters were found.  note that if (evtClass & SubScriberBlob.mask)
   * is true, it is considered to match.
   */
  bool unsubscribe(int evtClass, std::string streamName);

  /**
   * The most general logging class.  LogEvents are generic containers.  The
   * subscriber is given the event, and has to decide what to tod with it.
   * The is the most general case.  Several shortcuts are below.
   */
  inline void logit(int evtClass, const std::string& streamName, const LogEvent& e);

  /**
   * Log a simple message.
   */
  inline void logit(int evtClass, const std::string& streamName, const std::string& message);

  /**
   * Specialisations for the various event types.  They just call the above
   * two logit functions with appropriate presets.
   *
   * @{
   */
  inline void running (const std::string& streamName, const LogEvent& e);
  inline void moreInfo(const std::string& streamName, const LogEvent& e);
  inline void debug1  (const std::string& streamName, const LogEvent& e);
  inline void debug2  (const std::string& streamName, const LogEvent& e);
  inline void error   (const std::string& streamName, const LogEvent& e);
  inline void warn    (const std::string& streamName, const LogEvent& e);

  inline void running (const std::string& streamName, const std::string& message);
  inline void moreInfo(const std::string& streamName, const std::string& message);
  inline void debug1  (const std::string& streamName, const std::string& message);
  inline void debug2  (const std::string& streamName, const std::string& message);
  inline void error   (const std::string& streamName, const std::string& message);
  inline void warn    (const std::string& streamName, const std::string& message);
  /**
   * @}
   */

  /**
   * Check to see if a logging event will happen.  If logging an event
   * involves some work, it is worthwhile to check first, for efficiency's
   * sake.
   */
  bool willLogit(int evtClass, const std::string& streamName);

  /**
   * Specialisations for the various event types.  They just call the above
   * two logit functions with appropriate presets.
   *
   * @{
   */
  inline bool willRunning (const std::string& streamName);
  inline bool willMoreInfo(const std::string& streamName);
  inline bool willDebug1  (const std::string& streamName);
  inline bool willDebug2  (const std::string& streamName);
  inline bool willError   (const std::string& streamName);
  inline bool willWarn    (const std::string& streamName);
  /**
   * @}
   */

private:
  struct SubscriberBlob {
    ///
    std::string name;

    /// bitmap for which event types get logged.
    int mask;

    /// The actual subscriber
    RCPtr<LogWriter> sub;
  };

  /**
   * The list of subscribers.  We don't expect there to be more than 2 or
   * maybe 3. subscriptions don't happend often.
   */
  std::vector<SubscriberBlob> subscribers;

  /**
   * a bitmap that contains a summary of everything in the subscribers
   * list.  it is made up |ing all the subscribers' masks together to
   * enabling a very quick decision if nothing will be logged.
   */
  int mask;
};

// ============================= Globals ======================================

/**
 * This is a global logger for patrec to use.  The str() function in
 * stringops is very useful for calls to this.
 */
extern PRLITE_API Logger& getPatrecLogger();
#define PRLITE_pLog getPatrecLogger()

/** Subscribe a LogWriter, writing to the given filename, at the given
    level.

    These environment variables control it behaviour:
    LOGFILENAME
    LOGNAMES
    LOGLEVEL

  LOGFILENAME:

    Set filename to "-" for stderr instead of a file.  This is the default.
    Set it to a filename for redirecting the log output to a file.
    The logging system always appends the logfile.


  LOGLEVEL:

    Set it to 0,1,2,3 or 4. See the definition of Logger::LEVEL1... for more.
    Set it to 0 to switch logging off completely.
    The default is 1, which causes all output from LOGRUNNING, LOGWARN
    and LOGERROR to reach the output


  LOGNAMES:

    Subscribe only to specific names.  This must be a comma seperated list.

    For example: if you are only interested in what MyClass is doing,
    and MyClass logs using

      LOGRUNNING(pLog, isA(), "message");

    Then setting LOGNAMES=MyClass will do the trick.

    Examples:

    # to log to a file:
    export LOGFILENAME="/tmp/mylogfile"

    # to see all output put to LOGDEBUG2:
    export LOGLEVEL=4

    # to see nothing:
    export LOGLEVEL=0

    # to see output from these classes only:
    export LOGNAMES="CART_Tree,CART_Trainer,CART_QuestionFactory"

    Another way to set the variable for your program is to run it like
    so:

    LOGFILENAME="objfun.log" objfun

    and if you want to pipe your script through it:

    echo 0 | LOGFILENAME="objfun.log" objfun

*/
PRLITE_API void initLogging(const std::string& procname="process",
                            const std::string& filename="-");

// ============================= Inlines ======================================

void Logger::logit(int evtClass, const std::string& streamName, const LogEvent& e) {
  // quickquick test
  if (mask | evtClass) {
    // ok, test all
    std::vector<SubscriberBlob>::iterator iter = subscribers.begin();
    for ( ; iter!=subscribers.end(); iter++) {
      if ( (iter->mask & evtClass) &&
     ( (iter->name.find(streamName) == 0) || (iter->name == "") ) ) {
  iter->sub->logit(evtClass, streamName, e);
      }
    }
  }
}

void Logger::logit(int evtClass, const std::string& streamName, const std::string& message) {
  logit(evtClass, streamName, LogEvent(message));
}

void Logger::running (const std::string& streamName, const std::string& message) {
  logit(RUNNING_LOG, streamName, LogEvent(message));
}

void Logger::moreInfo(const std::string& streamName, const std::string& message) {
  logit(MORE_INFO, streamName, LogEvent(message));
}

void Logger::debug1  (const std::string& streamName, const std::string& message) {
  logit(DEBUG_LEVEL_1, streamName, LogEvent(message));
}

void Logger::debug2  (const std::string& streamName, const std::string& message) {
  logit(DEBUG_LEVEL_2, streamName, LogEvent(message));
}

void Logger::error   (const std::string& streamName, const std::string& message) {
  logit(ERROR_LOG, streamName, LogEvent(message));
}

void Logger::warn    (const std::string& streamName, const std::string& message) {
  logit(WARN, streamName, LogEvent(message));
}

void Logger::running (const std::string& streamName, const LogEvent& e) {
  logit(RUNNING_LOG, streamName, e);
}

void Logger::moreInfo(const std::string& streamName, const LogEvent& e) {
  logit(MORE_INFO, streamName, e);
}

void Logger::debug1  (const std::string& streamName, const LogEvent& e) {
  logit(DEBUG_LEVEL_1, streamName, e);
}

void Logger::debug2  (const std::string& streamName, const LogEvent& e) {
  logit(DEBUG_LEVEL_2, streamName, e);
}

void Logger::error   (const std::string& streamName, const LogEvent& e) {
  logit(ERROR_LOG, streamName, e);
}

void Logger::warn    (const std::string& streamName, const LogEvent& e) {
  logit(WARN, streamName, e);
}

bool Logger::willRunning (const std::string& streamName) {
  if (mask & RUNNING_LOG)
    return willLogit(RUNNING_LOG, streamName);
  else
    return false;
}

bool Logger::willMoreInfo(const std::string& streamName) {
  if (mask & MORE_INFO)
    return willLogit(MORE_INFO, streamName);
  else
    return false;
}

bool Logger::willDebug1  (const std::string& streamName) {
  if (mask & DEBUG_LEVEL_1)
    return willLogit(DEBUG_LEVEL_1, streamName);
  else
    return false;
}

bool Logger::willDebug2  (const std::string& streamName) {
  if (mask & DEBUG_LEVEL_2)
    return willLogit(DEBUG_LEVEL_2, streamName);
  else
    return false;
}

bool Logger::willError   (const std::string& streamName) {
  if (mask & ERROR_LOG)
    return willLogit(ERROR_LOG, streamName);
  else
    return false;
}

bool Logger::willWarn    (const std::string& streamName) {
  if (mask & WARN)
    return willLogit(WARN, streamName);
  else
    return false;
}

/*
 * Cool macros.  use these instead of calling the function directly.
 */
#define PRLITE_LOGDEBUG1(L,N,X) \
{ \
    if (L.willDebug1(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X; \
        L.debug1(N,oss.str()); \
    } \
}

#define PRLITE_LOGDEBUG2(L,N,X) \
{ \
    if (L.willDebug2(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X; \
        L.debug2(N,oss.str()); \
    } \
}

#define PRLITE_LOGRUNNING(L,N,X) \
{ \
    if (L.willRunning(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X; \
        L.running(N,oss.str()); \
    } \
}

#define PRLITE_LOGMOREINFO(L,N,X) \
{ \
    if (L.willMoreInfo(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X;  \
        L.moreInfo(N,oss.str()); \
    } \
}

#define PRLITE_LOGERROR(L,N,X) \
{ \
    if (L.willError(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X; \
        L.error(N,oss.str()); \
    } \
}

#define PRLITE_LOGWARN(L,N,X) \
{ \
    if (L.willWarn(N)) { \
        std::ostringstream oss; \
        oss << __FILE__ << ":" << __LINE__ << "|" << X; \
        L.warn(N,oss.str()); \
    } \
}

#define PRLITE_CAUTION(A)                                      \
  {                                                     \
    if (getPatrecLogger().willWarn("")) {               \
      std::ostringstream oss;                           \
      oss  << __FILE__ << ":" << __LINE__ << "|" << A;  \
      getPatrecLogger().warn("",oss.str());             \
    }                                                   \
  }


//} // namespace prlite

#endif // PRLITE_LOGGING_HPP
