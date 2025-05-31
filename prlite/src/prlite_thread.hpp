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

#ifndef PRLITE_THREAD_HPP
#define PRLITE_THREAD_HPP

/*
 * Author: Albert Visagie
 * Version: $Id$
 *
 * Purpose: Provide the very simplest synchronization primitives that
 * will allow building patrec in a thread safe way.  This is identical
 * to Phase 1 in the BOOST threading library.  As such it does not
 * provide methods to create threads.  That belongs outside patrec.
 *
 * We provide several implementations: one based on BOOST, one on
 * pthreads, one for Win32 (NIY) and a dummy one.  See thread.test.cc
 *
 * It provides only Mutex and RecursiveMutex, and locking
 * semantics. The idea is only to safe-guard things that will be
 * accessed by multiple threads, and not everything that might
 * possibly be.  Take for example HMMs: don't protect them, that
 * belongs at the application level.  We should however protect access
 * to the exemplar globals and the global logger.
 */


// patrec headers
#include "prlite_gendefs.hpp"  // PRLITE_API
#include "prlite_error.hpp"  // DEFINE_EXCEPTION

 //namespace prlite{

/**
 * Thrown for abuse of locking idioms.
 */
PRLITE_DEFINE_EXCEPTION(ThreadError, PatRecError);
PRLITE_DEFINE_EXCEPTION(MutexError, ThreadError);

//} // namespace prlite

#if HAVE_BOOST
#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/prlite_thread.hpp>

typedef boost::mutex Mutex;
typedef boost::recursive_mutex RecursiveMutex;

#elif THREAD_PTHREAD

#include <pthread.h>
#include <cstdlib>

//namespace prlite {

class PRLITE_API MutexBase {
protected:
  pthread_mutex_t my_mutex;
  pthread_mutexattr_t my_attr;
public:
  virtual ~MutexBase();
  virtual void lock();
  virtual void unlock();
};

class PRLITE_API ScopedLock {
  MutexBase* m;
public:
  ScopedLock(MutexBase& _m);
  ~ScopedLock();
};

class PRLITE_API Mutex : public MutexBase {
private:
  Mutex(const Mutex& e) : MutexBase(e) { }

public:
  Mutex();
  virtual void unlock();

  // this breaks the naming convention to remain boost compatible.
  typedef ScopedLock scoped_lock;
};


class PRLITE_API RecursiveMutex : public MutexBase {
private:
  RecursiveMutex(const RecursiveMutex& e) : MutexBase(e) { }

public:
  RecursiveMutex();

  virtual void unlock();

  // this breaks the naming convention to remain boost compatible.
  typedef ScopedLock scoped_lock;

};


#elif THREAD_WIN32

#error Win32 threads are not yet implemented.  Your best option is to use Boost.

#else // dummy threads

class PRLITE_API MutexBase {
public:
  void lock() { }
  void unlock() { }
  ~MutexBase() { }
};


class PRLITE_API ScopedLock {
public:
  ScopedLock(MutexBase&) { }
  ~ScopedLock() { }
};

class PRLITE_API Mutex : public MutexBase {
private:
  Mutex(const Mutex&) { }
public:
  Mutex() { }
  // this breaks the naming convention to make it boost compatible.
  typedef ScopedLock scoped_lock;
};

typedef Mutex RecursiveMutex;

//} // namespace prlite

#endif

#endif // PRLITE_THREAD_HPP
