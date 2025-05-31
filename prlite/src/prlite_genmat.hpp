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
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */


#ifndef PRLITE_GENMAT_HPP
#define PRLITE_GENMAT_HPP

// patrec headers
#include "prlite_genvec.hpp"  // Vector

// glinear headers
#include "gLinear/Matrix.h"  // gMatrix

// standard headers
#include <iosfwd>  // istream, ostream
#include <string>  // string

namespace prlite{

  template<typename T>
  using RowMatrix = gLinear::gMatrix<T, gLinear::ROW_DENSE>;

  template<typename T>
  using ColMatrix = gLinear::gMatrix<T, gLinear::COL_DENSE>;

  template<typename T>
  int length(const gLinear::gRowMatrix<T>& mat);

  template<typename T>
  class DiagMatrix;

/**
   ?
*/

  template<typename T>
  class DiagMatrix {

    //============================ Traits and Typedefs ============================

    //======================== Constructors and Destructor ========================
  public:

    DiagMatrix();
    DiagMatrix(const DiagMatrix<T>& m);
    DiagMatrix(int dim);
    explicit DiagMatrix(const gLinear::gRowMatrix<T>& m);
    ~DiagMatrix();

    //========================= Operators and Conversions =========================
  public:

    DiagMatrix<T>& operator=(const DiagMatrix<T>& m);
    DiagMatrix<T>& operator=(const T& elem);
    T& operator[](int i);
    T operator[](int i) const;
    operator gLinear::gRowMatrix<T>() const;

    // avoid these
    T operator()(int i) const;
    T& operator()(int i);
    T operator()(int r,
                 int c) const;
    T& operator()(int r,
                  int c);

    //================================= Iterators =================================
  public:

    typename gLinear::gRowVector<T>::T_ConstIterator begin() const;
    typename gLinear::gRowVector<T>::T_ConstIterator end() const;
    typename gLinear::gRowVector<T>::T_Iterator begin();
    typename gLinear::gRowVector<T>::T_Iterator end();

    //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

    //=========================== Inplace Configuration ===========================
  public:

    std::istream& read(std::istream& file);
    std::ostream& write(std::ostream& file) const;

    //===================== Required Virtual Member Functions =====================

    //====================== Other Virtual Member Functions =======================

    //======================= Non-virtual Member Functions ========================
  public:

    const std::string& isA() const;
    void resize(int sz);
    void deepen();
    void setAll(const T& val);
    T elem(int i) const;
    T& elem(int i);
    T elem(int r,
           int c) const;
    T& elem(int r,
            int c);
    int size() const;
    const T* arrayAddress() const;
    T* arrayAddress();
    void print(int depth = 0) const;

    //================================== Friends ==================================
  public:

    template<typename M>
    friend std::ostream& operator<<(std::ostream& file,
                                    const DiagMatrix<M>& dm);
    template<typename M>
    friend std::istream& operator>>(std::istream& file,
                                    DiagMatrix<M>& dm);

    //=============================== Data Members ================================
  private:

    gLinear::gRowVector<T> mainDiag;
    static T offDiag;    // off-diagonal values

    //============================ DEPRECATED Members =============================
  private:

    // T* odPtr;  // off-diagonal values

  }; // DiagMatrix


/**
   ?
*/

  template<typename T>
  class Transpose {

    //============================ Traits and Typedefs ============================

    //======================== Constructors and Destructor ========================
  public:

    Transpose(const gLinear::gRowMatrix<T>& mat);

    //========================= Operators and Conversions =========================
  public:

    T operator()(int i,
                 int j) const;
    T& operator()(int i,
                  int j);
    operator const gLinear::gRowMatrix<T>() const;

    //================================= Iterators =================================

    //============================= Exemplar Support ==============================
// there is nothing here because this class is not decended from ExemplarObject

    //=========================== Inplace Configuration ===========================

    //===================== Required Virtual Member Functions =====================

    //====================== Other Virtual Member Functions =======================

    //======================= Non-virtual Member Functions ========================
  public:

    int rows() const;
    int cols() const;
    T elem(int i,
           int j) const;
    T& elem(int i,
            int j);
    const std::string& isA() const;

    //================================== Friends ==================================

    //=============================== Data Members ================================
  private:

    gLinear::gRowMatrix<T> orig;

    //============================ DEPRECATED Members =============================
  public:

    // Vector<T>* arrayAddress() const;

  }; // Transpose

} // prlite

#include "prlite_genmat.tcc"



#endif // PRLITE_GENMAT_HPP
