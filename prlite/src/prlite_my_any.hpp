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
  $Id: prlite_my_any.hpp 4852 2006-10-23 20:02:03Z schwardt $
*/

/*
 * Author     : A.S. Visagie (DSP Group, E&E Eng, US)
 * Created on : 24 January 2005
 * Copyright  : University of Stellenbosch, all rights retained
 */

#ifndef PRLITE_MY_ANY_HPP
#define PRLITE_MY_ANY_HPP

// standard headers
#include <typeinfo>
#include <algorithm> // swap

// patrec headers
#include "prlite_error.hpp"

PRLITE_DEFINE_EXCEPTION(AnyBadCast, PatRecError);

//namespace prlite{

class Any;
std::ostream& operator<<(std::ostream& os, const Any& any);

/**
 * You can assign anything to it, and get it back. Getting it back is
 * dynamically type-checked.  If dynamic type checking fails, a
 * Any::exception is raised.  Any's can also be converted (by value!) to
 * their contents.  This is of course also type checked at run-time.
 *
 * Requirements on contained types:
 * 1.  They must support copy construction.
 * 2.  Copy construction must give a strong nothrow guarantee.
 *
 * Assigning of the content type is not required.
 *
 * see test_any.cc for example usage.
 *
 * Also responds to ABORT_NOT_THROW like gLinear and prlite_error.hpp.
 */
class Any {
private:
  // forwards
  class __HolderBase;
  template<typename T>
       class __Holder;

private:
  __HolderBase* contents;

public:

public:
  /// default constructor.
  Any() : contents(0) { }

  /// construct from contained type.
  template<typename T>
       Any(const T& t) : contents(new __Holder<T>(t)) { }

  /// copy constructor.  evetually uses the contained type's copy constructor to obtain our own copy of the object.
  Any(const Any& a) : contents(a.copy_contents()) { }

  /// destructor.  not virtual!
  ~Any() {
    if (contents) {
      delete contents;
    }
  }

  /// test for emptyness.
  bool empty(void) const {
    return contents == 0;
  }

  /**
   * print to a stream
   */
  std::ostream& print(std::ostream& os) const {
    return contents->print(os);
  }

  /**
   * get typeid( contents->contents ).  Can be used as
   *
   * Any any(123);
   * if (any.type() == typeid(int) { cout << "it's an int!" << endl; }
   *
   * if the above test succeeds, a cast<int>(any) will certainly not throw
   * exceptions, I think.
   */
  const std::type_info & type() const {
    return contents ? contents->type() : typeid(void);
  }

  /**
   * assign to arbitrary type.  Uses the type's copy constructor to get a copy.
   */
  template<typename T>
       Any& operator=(const T& t) {
   if (contents) {
     delete contents;
     contents = 0;
   }
   contents = new __Holder<T>(t);
   return *this;
       } // operator=

  /**
   * Assign to Any. self-assignment and exception safe.
   */
  Any& operator=(const Any& a) {
    Any tmp( a ); // make a copy.  may throw
    swap( tmp );  // our old contents go down with tmp.  can't throw.
    return *this;
  } // operator=

  /**
   * equalty.
   */
  bool operator==(const Any& a) const {
    if (contents && a.contents) {
      return contents->equal(a.contents);
    }
    // check for both being zero
    return ! (contents || a.contents);
  }

  /**
   * inequality.
   */
  bool operator!=(const Any& a) const {
    return ! operator==(a);
  }

  /**
   * Casting to contained type.  Does type checking. throws Any::exception in case.
   */
  template<typename T>
       T* get_contents(const T*) {
   if (!contents) {
     return 0;
   }

   __Holder<T>* g = dynamic_cast<__Holder<T>*>(contents);
   if ( g != NULL ) {
     return &g->get();
   }
   else {
     return 0;
   }
       }

  /**
   * equivalent(ish) to dynamic_cast<T*>(Any)
   */
  template<typename T>
       const T* get_contents(const T*) const {
   if (!contents) {
     return 0;
   }

   const __Holder<T>* g = dynamic_cast<const __Holder<T>*>(contents);
   if ( g != NULL ) {
     return &g->get();
   }
   else {
     return 0;
   }
       }

private:

  /**
   * Swap contents.
   */
  Any& swap(Any& other) throw() {
    std::swap(other.contents, contents);
    return *this;
  }


  /**
   * Helper to make copy constructor cleaner.
   */
  __HolderBase* copy_contents() const {
    if (contents) {
      return contents->copy();
    }
    else {
      return 0;
    }
  }

private:
  /**
   * Root class for the polymorphic container.
   */
  class __HolderBase {
  public:
    __HolderBase() { }
    virtual ~__HolderBase() { }

    // virtual copy
    virtual __HolderBase* copy() const = 0;

    // virtual equality
    virtual bool equal(const __HolderBase* b) const = 0;

    // type check
    virtual const std::type_info & type() const = 0;

    /**
     * print to a stream.
     */
    virtual std::ostream& print(std::ostream& os) const = 0;
  };

  /**
   * derived class that is instantiated to contain all the other types.
   * Implements some handles to deal with arbitrary types.
   */
  template<typename T>
       class __Holder : public __HolderBase
       {
       private:
         T contents;
       public:
         __Holder() { }
         __Holder(const T& t) : contents(t) { }
         __Holder(const __Holder& t) : __HolderBase(t), contents(t.contents) { }
         virtual ~__Holder() { }

         // to use these, the object must be dynamic_cast to __Holder<T> first.
         T& get() {
           return contents;
         }
         const T& get() const {
           return contents;
         }

         // virtual copy
         virtual __Holder* copy() const {
           return new __Holder(*this);
         }

         // virtual equality
         virtual bool equal(const __HolderBase* b) const {
           const __Holder* g;
           g = dynamic_cast<const __Holder*>(b);
           if (!g) {
             return false;
           }
           else {
             return g->contents == contents;
           }
         }

         /**
          * print to a stream.
          */
         virtual std::ostream& print(std::ostream& os) const {
           // Any'able types must have operator<<
           return os << contents;
         }

         // quicker than try { cast } type check
         virtual const std::type_info & type() const {
           return typeid(contents);
         }
       };
};


/**
 * utility function to get your type back out.  throws AnyBadCast if it
 * cannot cast.
 */
template<typename T>
     T& cast(Any& a) {
       T* t = NULL;
       t = a.get_contents(t);
       if (!t) {
   PRLITE_THROW(AnyBadCast, "");
       }
       return *t;
     }

/**
 * utility function to get your type back out.  throws AnyBadCast if it
 * cannot cast.
 */
template<typename T>
     const T& cast(const Any& a) {
  const T* t = NULL; // to shut the compiler up
  t = a.get_contents(t);
  if (!t) {
    PRLITE_THROW(AnyBadCast, "");
  }
  return *t;
}

/**
 * casting to a pointer to the internal type.
 */
template<typename T>
     T* pcast(Any& a) {
         T* t = NULL;
       return a.get_contents(t);
     }

/**
 * casting to a pointer to the contained type.
 */
template<typename T>
     const T* pcast(const Any& a) {
         T* t = NULL;
       return a.get_contents(t);
     }

/**
 * stream an Any out.
 */
inline std::ostream& operator<<(std::ostream& os, const Any& any) {
  return any.print(os);
}

/**
 * useful little helper: constructs an Any according to the type
 * string.
 */
Any mkAny(const std::string& value, const std::string& type);

/**
 * Another useful little helper: returns a type string depending on
 * the type of Any.
 */
std::string anyType(const Any& any);

//} // namespace prlite

#endif // PRLITE_MY_ANY_HPP
