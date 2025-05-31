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

char THREAD_CC_VERSION[] = "$Id$";

// patrec headers
#include "prlite_thread.hpp"
#include "prlite_error.hpp"

using namespace std;

#if THREAD_PTHREAD

MutexBase::~MutexBase() {
  int res = pthread_mutex_destroy(&my_mutex);
  pthread_mutexattr_destroy(&my_attr);
  if (res) {
    PRLITE_THROW(ThreadError, "pthread_mutex_destroy: mutex is still busy!");
  }
}

void MutexBase::lock() {
  //cerr << "lock() . . . 1" << endl;
  pthread_mutex_lock(&my_mutex);
  //cerr << "lock() . . . 2" << endl;
}

void MutexBase::unlock() {
  //cerr << "*** base called! ***" << endl;
}

ScopedLock::ScopedLock(MutexBase& _m) {
  //cerr << "locking . . . 1" << endl;
  _m.lock();
  //cerr << "locking . . . 2" << endl;
  m = &_m;
  //cerr << "locking . . . 3" << endl;
}

ScopedLock::~ScopedLock() {
  //cerr << "unlocking . . . 1" << endl;
  m->unlock();
  //cerr << "unlocking . . . 2" << endl;
}

Mutex::Mutex() {
  pthread_mutexattr_init(&my_attr);
  pthread_mutexattr_settype(&my_attr, PTHREAD_MUTEX_ERRORCHECK);
  pthread_mutex_init(&my_mutex, NULL);
}

void Mutex::unlock() {
  //cerr << "unlock() . . . 1" << endl;
  int res = pthread_mutex_unlock(&my_mutex);
  if (res != 0) {
    PRLITE_THROW(MutexError,
          "mutex unlocking error: it was either not locked, or not locked by the current thread.");
  }
  //cerr << "unlock() . . . 2" << endl;
}

RecursiveMutex::RecursiveMutex() {
  pthread_mutexattr_init(&my_attr);
  pthread_mutexattr_settype(&my_attr, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&my_mutex, &my_attr);
}

void RecursiveMutex::unlock() {
  //cerr << "unlock() . . . 1" << endl;
  pthread_mutex_unlock(&my_mutex);
  //cerr << "unlock() . . . 2" << endl;
}

#elif THREAD_WIN32

#error Win32 threads are not yet implemented.  Your best option is to use Boost.

#else

// all inline

#endif
