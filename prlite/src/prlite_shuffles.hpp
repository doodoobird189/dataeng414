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

#ifndef SHUFFLES_HPP
#define SHUFFLES_HPP

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14/09/95
         PURPOSE: Random sequence generator.

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"  // Vector
#include "prlite_mrandom.hpp"  // MRandom


namespace prlite{

/**
  The class ShuffledSequence supplies a facility for obtaining the indexes
  of a complete random sampling without replacement.
*/

class ShuffledSequence {

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
public:

  ShuffledSequence(int len);                         // len is no of indexes
  ShuffledSequence(int len,
                   long seed);

  //========================= Operators and Conversions =========================
public:

  int& operator[](int i);                            // index i

  //================================= Iterators =================================

  //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

  //=========================== Inplace Configuration ===========================

  //===================== Required Virtual Member Functions =====================

  //====================== Other Virtual Member Functions =======================

  //======================= Non-virtual Member Functions ========================
public:

  void randomize();
  void shuffle();                                    // new sampling order
  void resize(int len);                              // chg sample size
  void draw(gLinear::gRowVector<int>& dseq);                      // return next N from seq
  gLinear::gRowVector<int> draw(int len);
  int size() const;                                        // sample size
  void print(int depth = 0) const;                         // print indexes

  //================================== Friends ==================================

  //=============================== Data Members ================================
private:

  gLinear::gRowVector<int> seq;
  int alreadyDrawn;
  prlite::MRandom gen;

  //============================ DEPRECATED Members =============================


}; // ShuffledSequence

} // prlite

#endif // SHUFFLES_HPP
