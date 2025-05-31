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

char SHUFFLES_CC_VERSION[] = "$Id$";

/*******************************************************************************

         AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
         DATE:    14/09/95
         PURPOSE: Random sequence generator.

*******************************************************************************/

// patrec headers
#include "prlite_shuffles.hpp"
#include "prlite_error.hpp"

// standard headers
#include <iomanip>

using namespace std;

namespace prlite{

/******************************************************************************/

ShuffledSequence::ShuffledSequence(int len)
  : seq(len) {
  shuffle();
  alreadyDrawn = 0;
} // ShuffledSequence

ShuffledSequence::ShuffledSequence(int len,
                                   long seed)
  : seq(len),
    gen(seed) {
  shuffle();
  alreadyDrawn = 0;
} // ShuffledSequence

void ShuffledSequence::randomize() {
  gen.randomize();
} // randomize

void ShuffledSequence::shuffle() {
  gLinear::gRowVector<int> available( seq.size() );
  int i;
  for (i = 0; i < available.size(); i++) {
    available[i] = i;
  } // for
  int cnt = seq.size();
  for (i = 0; i < seq.size(); i++) {
    int choose = gen.Int(cnt);
    seq[i] = available[choose];
    available[choose] = available[cnt - 1];
    cnt--;
  } // for i
  alreadyDrawn = 0;
} // shuffle

void ShuffledSequence::resize(int len) { // chg sample size
  seq.resize(len);
  shuffle();
} // resize

int& ShuffledSequence::operator[](int i) {           // index i
  return seq[i];
} // operator[]

void ShuffledSequence::draw(gLinear::gRowVector<int>& dseq) {
  if ( alreadyDrawn + dseq.size() > seq.size() ) {
    PRLITE_FAIL("ShuffledSequence::draw : cannot draw " << dseq.size()
         << " values, only " << (seq.size() - alreadyDrawn) << " available.");
  } // if
  for (int i = 0; i < dseq.size(); i++) {
    dseq[i] = seq[i + alreadyDrawn];
  } // for i
  alreadyDrawn += dseq.size();
} // draw

gLinear::gRowVector<int> ShuffledSequence::draw(int len) {
  if ( alreadyDrawn + len > seq.size() ) {
    PRLITE_FAIL("ShuffledSequence::draw : cannot draw " << len
         << " values, only " << (seq.size() - alreadyDrawn) << " available.");
  } // if
  gLinear::gRowVector<int> dseq(len);
  for (int i = 0; i < dseq.size(); i++) {
    dseq[i] = seq[i + alreadyDrawn];
  } // for i
  alreadyDrawn += dseq.size();
  return dseq;
} // draw

int ShuffledSequence::size() const {
  return seq.size();
} // size

void ShuffledSequence::print(int depth) const {
  if (depth >= 0) {
    for (int i = 0; i < seq.size(); i++) {
      cout << setw(8) << seq[i];
    } // for
    cout << endl;
  } // if
} // print

} // prlite
