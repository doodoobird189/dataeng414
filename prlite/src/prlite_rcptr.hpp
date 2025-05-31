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

/*******************************************************************************
AUTHOR: J.A du Preez (Copyright University of Stellenbosch, all rights reserved.)
*******************************************************************************/

#ifndef PRLITE_RCPTR_HPP
#define PRLITE_RCPTR_HPP

// standard headers
#include <iosfwd>  // ostream
#include <ciso646>  // and etc

 //namespace prlite{

/**
This class implements a Reference Counted Pointer to an object -
it is used similar to a built-in pointer, but is much safer.
Semantically it is related to the "handle" used in JAVA.\\

{\bfADVANTAGES over standard pointers:}\\

(a) Copying or assigning RCPtr's provide safe control over the
lifetime of the commonly pointed to object. (Only the last surviving
RCPtr will destroy the common object.) No multiple deletes or
who's the boss problems!\\

(b) You never explicitly delete memory pointed to by this class, it will
automatically apply a type of garbage collection. You use it and forget
it. No memory leaks!

(c) RCPtr's are initialised to 0 if nothing else is provided.
No dangling pointers!

(d) The ++ and -- operations are disabled, this is a pointer, not an array!\\

{\bfDIFFERENCES with standard pointers:}\\

(a) Conversion to higher level (less derived) RCPtr's are done via constructors
which creates temporary objects pointing to the same data. This is typically transparent
but consider:\\

    void func( RCPtr<Base>& ) {implementation...}\\

    RCPtr<Derived> d; func(d) \\

Problem is that a temporary RCPtr<Base> will be created. Although this is fine,
the compiler will generate a warning. If you want to avoid it, first construct/assign to
a RCPtr<Base>. (No problems with const references).\\

(b) Conversion to lower level (more derived) or conversion from const to non-const RCPtr's
needs to be done explicitly by using the standalone template functions 'dynamicCast', and/or
'constCast' (described elsewhere). NOTE - the small difference in the use of these functions:
dynamicCast<DestType>(sourcePtr) and NOT dynamicCast< RCPtr<DestType> > which is a closer parallel
to the behaviour of the standard dynamic_cast function and similarly for constCast.\\

(c) There is no support for covariant pointer returns - the requirement for it should be
fairly rare, covariant pointer returns are only essential when you use multiple
inheritance and more than one of the multiple parents specify the same virtual function
but with different pointer returns. Use standard pointers where this is essential.
The following will go totally awry when you return RCPtr's - the compiler will (correctly)
simply and utterly refuse to do this (if you don't know what this is, go and read it up or ask). \\

    class Base {
      virtual RCPtr<Base> func(parameters) = 0;
    }; // Base\\

    class Derived : public Base {
      virtual RCPtr<Derived> func(same parameters) {implementation}
    }; // Derived\\

It simply does not know of the implied inheritance relationship between the RCPtr's.
Either use RCPtr<Base> as returne type throughout, or return a naked Derived*
(controlling its own memory) in the derived class (which can then be swallowed by
another RCPtr shortly afterwards).\\

NOTE -  When the info you want returned is already contained in a RCPtr and you
want to transfer that to another RCPtr, detouring via a naked pointer pointing
to a (newly created) copy of the info is somewhat inefficient. (As already discussed
in Warning 1, if you're not working with a new copy and you're doing this, you are
in deep trouble in any case.)\\

WARNING 1: \\
Take the WARNINGS in the RCPtr member documentation (see in class) seriously,
if not your code will CRRROAK in the most horribly spectacular fashion!
In general, when you give a RCPtr control over a naked pointer, that RCPtr, together
with the other RCPtr's that it is AWARE of (i.e. they share a common ref count),
gets ALL control and responsibility for ever and ever amen. A sort of a Kremlin inner
ring of power. NOTE - RCPtr's pointing to the same data are not necessarily AWARE of
each other, that depends on how they were created/assigned. If they are not you will
regret it.) Woe be to you if something else destroys the naked pointer, or if the RCPtr
departs to code heaven too early. That also holds for abusing naked pointers obtained
via the objPtr() and cntPtr members. Be especially wary about the following types of
situations:\\

    {
      RCPtr<ClassType> aRCPtr;
      {
        ClassType anObject;
        aRCPtr = RCPtr<ClassType>(&anObject);
      }
      When you get here your 'aRCPtr' will be dangling in deep space - (of course
      this is no diffent from normal pointers)
    }\\

    or just as bad:\\

    {
      ClassType anObject;
      {
        RCPtr<ClassType> aRCPtr = RCPtr<ClassType>(&anObject);
      }
      When you get here your 'anObject' will be insane
    }\\
NOTE - In general using the address of (&) operator with RCPtr's is bad news.\\

WARNING 2: Avoid initialising temporary RCPtr's with naked pointers - your naked pointer will
get nuked. (Its playing with the big boys now.) When using RCPtr's as parameters to functions,
keep your eyes peeled for the following gotcha. Say you have a function:\\

    func(const RCPtr<ClassType>& par)\\

startling things may happen if you call it with a naked pointer ie\\

    {
      ClassType* nakedPtr = new ClassType;
      func( RCPtr<ClassType>(nakedPtr) );
    }\\

A temporary RCPtr is created when calling the function. It will destroy nakedPtr when it
auto-destructs, leaving nakedPtr dangling in deep space. NOTE - you have to give your
permission for this mayhem to happen. The constructor using a naked pointer has been
declared as "explicit" (see code). In other words, this constructor will not happen
without you explicitly saying so. And that is the ONLY way for a lonely naked pointer
to get into a RCPtr. Therefore you are fairly safe, it will not sneak up on you - quite
often it is the efficient option to use an existing RCPtr as a parameter to a function.\\

WARNING 3: In spite of all the nice things claimed of the RCPtr, you can still get a
memory leak even while using them! This happens when you get your RCPtrs tied up in a
loop (a doubly linked list using RCPtrs is one example). Here is the simplest example,
more elaborate ones are possible too (sketching it may help):\\

    class Type {
    public:
      RCPtr<Type> aPtr;
    }; // Type\\

Now (mis)use this class in a piece of code as follows:\\

    RCPtr<Type> extPtr(new Type);
    extPtr->aPtr = extPtr; // both pointers now look at the same object;
    extPtr = 0; // this one's a goner
    // but the internal aPtr won't destruct, it is still looking at something.
    // in other words there's now a memory leak here.\\

Another variant:\\

    RCPtr<Type> ext1Ptr(new Type);
    RCPtr<Type> ext2Ptr(new Type);
    ext1Ptr->aPtr = ext2Ptr;
    ext2Ptr->aPtr = ext1Ptr; // our pointer loop is now in place
    ext1Ptr = ext2Ptr = 0;   // gone with the wind
    // once again the two aPtr's won't destruct - memory leak!\\

How to fix it? Ok, this is tricky and situation dependend - you got to think carefully
about what you are doing. Have a good look at your loop of RCPtrs and see whether you
can break the loop by replacing certain of them with plain pointers (Yea, I know, you
thought you'll never use them again. Well ...) For instance, in a doubly linked list,
how about replacing the backward links with plain C pointers? That sort of thing.\\

TIP 1: You want to get a naked pointer out of the grips of a RCPtr? Currently there is
no direct support for this. And I can see no safe way of implementing it except to
tell the RCPtr's that they have just experienced a coup-de-tat and must go and sit and
sulk in the corner (or face a firing squad). These types of things, however, tend to
get gory and really mess existing states up. However, if you use the objPtr() member
you can get to the naked pointer, just remember that you may not control its lifetime.
And you can get a NEW/safe naked pointer by using a copy constructor to do
something like:\\

    ClassType* nakedPtr = new ClassType(*aRCPtr);\\

or if something like a copy/clone member is supported\\

    ClassType* nakedPtr = rcPtr->copy();\\

TIP 2:
As with all pointers, copies of any sort are shallow. In keeping with a minimalist
spirit, no support is provided for deep copies. To clone a newly constructed deep
object from an RCPtr, use the copy constructor of the underlying object on the dereference eg: \\
    T deep(*shallowRCPtr);           // construct from underlying object \\
Remember that if the copy constructor also was shallow you will have to deepen that too ie:
    deep.deepen();\\
You will often need to deepen an existing RCPtr. If it points to a standard patrec object that
sports a copy and deepen member, use the DEEPEN macro (Provided lower down in this file). With
non-standard patrec objects you can, similar to the above, use something like\\
RCPtr<T> deepPtr( RCPtr<T> ( new T(*shallowPtr) ) ); \\

TIP 3:
Play with rcptr.eg.cc for an example of typical use. It compiles standalone.

*/

template<typename T>
class RCPtr {

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
public:

  /// Default constructor, points to an nil object
  RCPtr();

  /**
   * Copy constructor, both RCPtr's now point to the same underlying object
   * and are aware of each other.
   */
  RCPtr(const RCPtr& rcp);

  /**
   * Copy constructor from a lower type RCPtr, both now point to the same underlying object
   * and are aware of each other.
   */
  template<typename FromType>
  RCPtr(const RCPtr<FromType>& rcp);

  /// Destructor, will not destroy the actual object if another RCPtr references it.
  ~RCPtr();

  /**
   * {\em WARNING:} Takes over lifetime control. You hereby solemnly swear that nothing external will
   * ever try to destroy justNewed. RCPtr may, however, destroy it whenever it deems necessary.
   * NOTE: Via this constructor the RCPtr has no knowledge of anything else (including other
   * RCPtr's) that may also lay claim to justNewed - a sure recipe for chaos! Can construct with NIL.
   */
  explicit RCPtr(T* justNewed);

  /**
   * Low-level constructor to Share the pointed-to object with another (compatible)
   * reference counting class.
   * The RCPtr and the other ref count class now point to the same underlying object
   * and are aware of each other. Useful for implementing a type of upcasting between
   * RCPtr's, but also allows ref counted sharing with potentially quite different
   * ref count classes. Recommended to prefer higher-level constructors.
   */
  RCPtr(T* objectPtr,
        unsigned* countPtr);

  //========================= Operators and Conversions =========================
public:

  /**
   * Assignment, both RCPtr's now point to the same underlying object and
   * are aware of each other.
   */
  RCPtr& operator=(const RCPtr& rhs);

  /**
   * Assignment, both RCPtr's now point to the same underlying object and
   * are aware of each other.
   */
  template<typename FromType>
  RCPtr<T>& operator=(const RCPtr<FromType>& rhs);

  /// Very handy to access all the members of the underlying object via ->
  T* operator->() const;

  /// Dereferencing also works
  T& operator*() const;

  /**
   * Less than operator, usefull for ordering RCPtrs. Comparison on
   * the underlying pointer itself, not the value of what is being
   * pointed to. Useful when including RCPtrs in associative data
   * structures
   */
  bool operator<(const RCPtr& rhs) const;

  /**
   * Equality operator. Comparison on the underlying pointer itself,
   * not the value that is being pointed to. Useful when including
   * RCPtrs in associative data structures
   */
  bool operator==(const RCPtr& rhs) const;

  /**
   * Inequality operator. Comparison on the underlying pointer itself,
   * not the value that is being pointed to. Useful when including
   * RCPtrs in associative data structures
   */
  bool operator!=(const RCPtr& rhs) const;

  /// For boolean comparisons, checks whether it has been instantiated
  operator bool() const;

  /**
   * Low-level assignment from another compatible reference counted object
   * (including RCPtr's), both now point to the same underlying object and
   * are aware of each other. Recommend to prefer the (safer) higher-level
   * members.
   */
  RCPtr& assignFrom(T* objectPtr,
                    unsigned* countPtr);

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

  //=========================== Inplace Configuration ===========================

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================

  //======================= Non-virtual Member Functions ========================
public:

  /**
   * Return the pointer that this object encapsulates.
   * {\em WARNING:} Please do not pass this returned pointer around
   * or abuse (e.g. delete) it. If you don't understand the issues
   * involved then rather refrain from using this call or ask someone
   * to explain it to you.
   */
  T* objPtr() const;

  /**
   * Return the counts pointer containing the reference counts.
   * {\em WARNING:} Please do not pass this returned pointer around
   * or abuse (e.g. delete) it. If you don't understand the issues
   * involved then rather refrain from using this call or ask someone
   * to explain it to you.
   */
  unsigned* cntPtr() const;

  /**
   * The pointed to object will now never be destroyed by the RCPtr. This
   * is useful if you want to assign the address of an existing object,
   * BUT {\em WARNING:} USE WITH EXTREME CARE, taking the address of an
   * existing object is dangerous stuff. The onus is now on YOU to make
   * absolutely sure that the existing object will not expire too early
   * (dangling RCPtr) or too late (memory leak)! Also note, this method
   * only holds for the currently pointed to object. Should you assign
   * a new one, it is a whole new ball-game.
   */
  void makeImmortal();

  /// true if object shared with other RCPtr's
  bool shared() const;

  //================================== Friends ==================================
public:

  /**
   * Dynamically checked cast to different type RCPtr.
   * Usage: toTypePtr = dynamicCast<ToType>(fromPtr)
   * NOTE - caller needs to check for success
   * @param fromPtr The source type RCPtr.
   * @return The casted RCPtr - 0 if failed.
   */
  template<typename ToType, typename FromType>
  friend RCPtr<ToType> dynamicCast(const RCPtr<FromType>& fromPtr);

  /**
   * Casting constness away from a RCPtr.
   * Usage: nonConstPtr = constCast<ToType>(fromPtr)
   * @param fromPtr The source type RCPtr.
   * @return The const casted RCPtr.
   */
  template<typename ToType, typename FromType>
  friend RCPtr<ToType> constCast(const RCPtr<FromType>& fromPtr);

  /**
   * Dynamically checked assignment implementing toPtr = fromPtr inplace.
   * Remember: assignTo(lhs, rhs) as in lhs = rhs;
   * NOTE - caller needs to check for success
   * Slightly more efficient than a dynamicCast operator= combination.
   * @param toPtr The destination type RCPtr.
   * @param fromPtr The source type RCPtr.
   * @return toPtr - The casted RCPtr - 0 if failed.
   * @return true if succesful
   */
  template<typename ToType, typename FromType>
  friend bool assignTo(RCPtr<ToType>& toPtr,
                       const RCPtr<FromType>& fromPtr);

  /// Can output it
  template<typename M>
  friend std::ostream& operator<<(std::ostream& file,
                                  const RCPtr<M>& rcp);

  //=============================== Data Members ================================
private:

  T* oPtr;
  unsigned* cPtr;

  //============================ DEPRECATED Members =============================


}; // RCPtr

#include "prlite_rcptr.tcc"

//} // namespace prlite

//namespace prlite{

#define PRLITE_DEEPEN(X,TYPE)                   \
if (X.shared()) {                               \
    X = RCPtr<TYPE>(X->copy());                 \
  }                                             \
  if (X) {                                      \
    X->deepen();                                \
  }

//} // namespace prlite

#endif // PRLITE_RCPTR_HPP
