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

char TESTNEW_CC_VERSION[] = "$Id$";

//#define TEST_SPECIAL

#include "gLinear.h"

#include <iostream>
using namespace std;
using namespace gLinear;

#define N 25
#define SRS 3
#define SRE 19
#define STR 3

#define LNR 8
#define LNC 10

#define NR 4
#define NC 5

#define MSRS 1
#define MSRE 6
#define MSRT 2
#define MSCS 2
#define MSCE 8
#define MSCT 3

// print separator
void separate(void) {
   for (int i = 0; i < 50; i++) {
      cout << "-";
   }
   cout << endl;
}

/* General set generating functions */



/* Vector */

// generic vector dump: using operator[]
template<typename T_Vector>
void vec_dump(const T_Vector& vec) {
  for (typename T_Vector::T_Index indx = 0; indx < vec.size(); indx++) {
    cout << vec[indx] << " ";
  }
  cout << endl;
}

// generic vector dump: using iterator
template<typename T_Vector>
void vec_iter_forw(const T_Vector& vec) {
  // create 'shortcut' type to the argument vector type iter
  const typename T_Vector::T_ConstIterator vecend = vec.end();
  for (typename T_Vector::T_ConstIterator iter = vec.begin(); iter != vecend; iter++) {
    cout << *iter << " ";
  }
  cout << endl;
}

// generic vector backward dump: using iteration
template<typename T_Vector>
void vec_iter_backw(const T_Vector& vec) {
  // create 'shortcut' type to the argument vector iterator type
  typedef typename T_Vector::T_ConstIterator T_Iterator;
  const T_Iterator vecbegin = vec.begin();
  for (T_Iterator iter = vec.end(); iter-- != vecbegin; ) {
    cout << *iter << " ";
  }
  cout << endl;
}

// generic vector set assign
template<typename T_Vector>
void vec_set(T_Vector vec, typename T_Vector::T_Value start = 0,
       typename T_Vector::T_Value step = 1) {
  typedef typename T_Vector::T_Value T_Value;
  typedef typename T_Vector::T_Iterator T_Iterator;
  const T_Iterator vecend = vec.end();
  T_Value i = start;
  for (T_Iterator iter = vec.begin(); iter != vecend; iter++) {
    *iter = i;
    i += step;
  }
}

/* Matrix */

// Edward dV 23/05/2002 -- changed return type from T_Matrix& to void to remove compiler warnings
// generic matrix set assign
template<typename T_Matrix>
void mat_set(T_Matrix mat,
       typename T_Matrix::T_Value start = 0,
       typename T_Matrix::T_Value step = 1) {
  typedef typename T_Matrix::T_Value T_Value;
  T_Value i = start;
  for (typename T_Matrix::T_RowIterator riter = mat.rowBegin(); riter != mat.rowEnd(); riter++) {
    vec_set(*riter, i, step); // use vector set assign
    i += static_cast<T_Value>( static_cast<float>( mat.cols() ) )  * step; //extra cast added to remove compiler warning. Not sure why it works (or complained).
  }
}

// generic matrix dump: using operator()
template<typename T_Value, typename T_Storage>
void mat_dump(const gMatrix<T_Value, T_Storage>& mat) {
  typedef gMatrix<T_Value, T_Storage> T_Matrix;
  // using indices
  for (typename T_Matrix::T_Index rindx = 0; rindx < mat.rows(); rindx++) {
    for (typename T_Matrix::T_Index cindx = 0; cindx < mat.cols(); cindx++) {
      cout << mat(rindx, cindx) << " ";
    }
    cout << endl;
  }
  cout << endl;
}

// generic matrix dump: using row iteration
template<typename T_Matrix>
void mat_iter_row(const T_Matrix& mat) {
  typedef typename T_Matrix::T_Row T_Row;
  typedef typename T_Matrix::T_ConstRowIterator T_RowIterator;
  // using row indexing and column iteration
  const T_RowIterator row_end = mat.rowEnd();
  for (T_RowIterator riter = mat.rowBegin(); riter != row_end; riter++) {
    vec_dump(*riter); // use vector dump
  }
}

// generic transpose matrix dump: using column iteration
template<typename T_Matrix>
void mat_iter_col(const T_Matrix& mat) {
  typedef typename T_Matrix::T_Col T_Col;
  typedef typename T_Matrix::T_ConstColIterator T_ColIterator;
  // using row indexing and column iteration
  const T_ColIterator col_end = mat.colEnd();
  for (T_ColIterator citer = mat.colBegin(); citer != col_end; citer++) {
    vec_dump(*citer); // use vector dump
  }
  }

// TMP to do sum(i=1 .. N) i
template<int N_Val>
struct sumOfN {
  enum {
    val = (N_Val + sumOfN<N_Val-1>::val)
  };
};
// spcialize on stop condition
template<>
struct sumOfN<1> {
  enum { val = 1
  };
};

// TMP to do sum(i=1..N) i^2
template<int N_Val>
struct sumOfNSq {
  enum {
    val = ( N_Val*N_Val + sumOfNSq<N_Val-1>::val )
  };
};
// stop condition
template<>
struct sumOfNSq<1> {
  enum {
    val = 1
  };
};

template<typename T_Vector>
void vec_testctor(const T_Vector&) {

  int r = 1;
  separate();
  cout << "ctors" << endl;
  separate();
  cout << "default ctor (no args)" << endl;
  T_Vector vec;
  cout << "vec.size(): " << vec.size() << " | " << 0 << endl;
  cout << endl;

  separate();
  cout << "dtor" << endl;
  vec.~T_Vector();
  cout << "ctor specifying size (using placement new)" << endl;
  // try out the placement new
  new(&vec) T_Vector(N);
  cout << "vec.size(): " << vec.size() << " | " << N << endl;
  cout << "initialize" << endl;
  vec_set(vec);
  cout << "vec:" << endl;
  vec_dump(vec);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "copy ctor: vec2(vec)" << endl;
  T_Vector vec2(vec);
  r++;
  cout << "vec2.size(): " << vec2.size() << " | " << N << endl;
  cout << "vec, vec2:" << endl;
  vec_dump(vec);
  vec_dump(vec2);
  cout << "vec2.getRefs(): " << vec2.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "creating and init'ing vec3 (offset 2)" << endl;
  T_Vector vec3(N);
  vec_set(vec3, 2);
  cout << "vec, vec3:" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "aliasing vec3: vec3.alias(vec): " << endl;
  vec3.alias(vec);
  r++;
  cout << "vec, vec3:" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vec3.getRefs(): " << vec3.getRefs() << " | " << r << endl;
  cout << "vec += 3:" << endl;
  vec += 3;
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "deepening vec3: vec3.deepen(): " << endl;
  vec3.deepen();
  r--;
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vec3.getRefs(): " << vec3.getRefs() << " | " << 1 << endl;
  vec -= 2;
  vec3 += 2;
  cout << "vec -= 2, vec3 += 2" << endl;
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);

  separate();
  cout << "assignment: vector-vector" << endl;
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vec3.getRefs(): " << vec3.getRefs() << " | " << 1 << endl;
  cout << "vec3 = vec" << endl;
  vec3 = vec;
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vec3.getRefs(): " << vec3.getRefs() << " | " << 1 << endl;
  cout << "assignment: vector-scalar" << endl;
  cout << "vec3 = 3" << endl;
  // JADP 17 Maart 2002 vec3 = 3;
  vec3.assignToAll(3);
  cout << "vec, vec3" << endl;
  vec_dump(vec);
  vec_dump(vec3);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vec3.getRefs(): " << vec3.getRefs() << " | " << 1 << endl;

  separate();
  cout << "vec: " << endl;
  vec_dump(vec);
  cout << "transpose: vectrans = vec.transpose()" << endl;
  typename T_Vector::T_Transpose vectrans = vec.transpose();
  r++;
  cout << "vectrans.size(): " << vectrans.size() << " | " << vec.size() << endl;
  cout << "vec, vectrans:" << endl;
  vec_dump(vec);
  vec_dump(vectrans);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << "vectrans.getRefs(): " << vectrans.getRefs() << " | " << r << endl;
}

template<typename T_Vector>
void vec_testsub(const T_Vector&) {

  int r = 1;
  typedef typename T_Vector::T_Slice T_Slice;
  typedef typename T_Vector::T_Subrange T_Subrange;
  typedef typename T_Vector::T_Index T_Index;

  T_Vector vec(N);
  vec_set(vec);

  separate();
  cout << "subranges and slices" << endl;
  separate();
  cout << "vec:" << endl;
  vec_dump(vec);
  cout << "iterating vector: " << endl;
  vec_iter_forw(vec);
  vec_iter_backw(vec);
  cout << "vec.getRefs(): " << vec.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "subranging vec: " << SRS << ":" << SRE << endl;
  T_Subrange vecsub = vec.subrange(SRS, SRE);
  r++;
  cout << "vecsub.size(): " << vecsub.size() << " | " << SRE - SRS + 1 << endl;
  cout << "vecsub:" << endl;
  vec_dump(vecsub);
  cout << "iterating vecsub:" << endl;
  vec_iter_forw(vecsub);
  vec_iter_backw(vecsub);
  cout << "vecsub.getRefs(): " << vecsub.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "slicing vec: " << SRS << ":" << STR << ":" << SRE << endl;
  T_Slice vecslice = vec.slice( gIndexRange(SRS, SRE, STR) );
  r++;
  cout << "vecslice.size(): " << vecslice.size() << " | " << ceil( static_cast<double>(SRE - SRS + 1) / static_cast<double>(STR) ) << endl;
  cout << "vecslice:" << endl;
  vec_dump(vecslice);
  cout << "iterating vecslice" << endl;
  vec_iter_forw(vecslice);
  vec_iter_backw(vecslice);
  cout << "vecslice.getRefs(): " << vecslice.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "subranging vecslice: " << "1..4" << endl;
  T_Slice vecsub2 = vecslice.subrange(1,4);
  r++;
  cout << "vecsub2.size(): " << vecsub2.size() << " | " << (4 - 1) + 1 << endl;
  cout << "vecsub2:" << endl;
  vec_dump(vecsub2);
  cout << "iterating vecsub2" << endl;
  vec_iter_forw(vecsub2);
  vec_iter_backw(vecsub2);
  cout << "vecsub2.getRefs(): " << vecsub.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "slicing vecslice: " << 1 << ":" << 2 << ":" << 4 << endl;
  T_Slice vecslice2 = vecslice.slice( gIndexRange(1,4,2) );
  r++;
  cout << "vecslice2.size(): " << vecslice2.size() << " | "
       << std::ceil( static_cast<float>((4 - 1 + 1) / 2) ) << endl;
  cout << "vecslice:" << endl;
  vec_dump(vecslice2);
  cout << "iterating vecsclice" << endl;
  vec_iter_forw(vecslice2);
  vec_iter_backw(vecslice2);
  cout << "vecslice2.getRefs(): " << vecslice2.getRefs() << " | " << r << endl;
  cout << endl;

  // abuse of aliasing
  cout << "nasty slicing" << endl;
  separate();
  cout << "slice = slice" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  vec.slice( gIndexRange(13,17,4) ) = vecslice2;
  cout << " 13,17 -> 6,12" << endl;
  cout << "vec: " << endl;
  vec_dump(vec);
  cout << endl;

  separate();
  cout << "subrange = slice" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  vec.subrange(2,7) = vecslice;
  cout << " 2..7 -> 3:3:18" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  cout << endl;

  separate();
  cout << "slice = subrange" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  vec.slice( gIndexRange(19,23,2) ) = vec.subrange(8,10);
  cout << " 19:2:23 -> 8..10" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  cout << endl;
}

template<typename T_Vector>
void vec_testmemopers(const T_Vector&) {

  T_Vector vec(N);
  T_Vector vec2(N);
  vec.assignToAll(0);
  vec_set(vec2);

  cout << "member operators & functions" << endl;
  separate();
  cout << "vec, vec2" << endl;
  vec_dump(vec);
  vec_dump(vec2);
  cout << "assignment: vector-vector" << endl;
  cout << "vec = vec2" << endl;
  vec = vec2;
  vec_dump(vec);
  vec_dump(vec2);

  separate();
  cout << "assignment: vector-scalar" << endl;
  // JADP 17 Maart 2002 vec = 3;
  vec.assignToAll(3);
  cout << "vec = 3:" << endl;
  vec_dump(vec);
  vec_set(vec);
  cout << "scalar computed (+ - * /) assignment" << endl;
  cout << "vec:" << endl;
  vec_dump(vec);
  vec += 1;
  cout << "vec += 1:" << endl;
  vec_dump(vec);
  vec -= 1;
  cout << "vec -= 1:" << endl;
  vec_dump(vec);
  vec *= 2;
  cout << "vec *= 2:" << endl;
  vec_dump(vec);
  vec /= 2;
  cout << "vec /= 2:" << endl;
  vec_dump(vec);
  cout << "vector computed (+ -) assignment" << endl;
  vec += 2;
  vec2 = 3*vec;
  cout << "vec, vec2:" << endl;
  vec_dump(vec);
  vec_dump(vec2);
  vec += vec2;
  cout << "vec += vec2:" << endl;
  vec_dump(vec);
  vec -= vec2;
  cout << "vec -= vec2:" << endl;
  vec_dump(vec);

  cout << "sum" << endl;
  vec_set(vec);
  vec += 1;
  cout << "vec: " << endl;
  vec_dump(vec);
  cout << "vec.sum(): " << vec.sum() << " | " << static_cast<int>(sumOfN<N>::val) << endl;
  cout << "vec.sum2(): " << vec.sum2() << " | " << static_cast<int>(sumOfNSq<N>::val) << endl;

}

template<typename T_Vector>
void vec_testopers(const T_Vector&) {

  T_Vector vec(N);
  T_Vector vec2(N);
  cout << "unary vector operators" << endl;
  separate();
  vec_set(vec);
  cout << "vec:, -vec:" << endl;
  vec_dump(vec);
  vec_dump(-vec);
  separate();
  cout << "global operators" << endl;
  separate();
  cout << "vec-vec" << endl;
  vec_set(vec);
  vec_set(vec2);
  cout << "equality and inequality" << endl;
  cout << "vec:, vec2: " << endl;
  vec_dump(vec);
  vec_dump(vec2);
  cout << (vec == vec2) << " " << (vec != vec2) << " | " << "1 0" << endl;
  vec_set(vec2, 3);
  cout << "vec:, vec2: " << endl;
  vec_dump(vec);
  vec_dump(vec2);
  cout << (vec == vec2) << " " << (vec != vec2) << " | " << "0 1" << endl;
  cout << "vec + vec2" << endl;
  vec_dump( (vec + vec2) );
  cout << "vec - vec2" << endl;
  vec_dump( (vec - vec2) );

  separate();
  cout << "vec-scalar" << endl;
  vec_set(vec);
  cout << "vec: " << endl;
  vec_dump(vec);
  cout << "vec + 2" << endl;
  vec_dump(vec + 2);
  cout << "vec - 2" << endl;
  vec_dump(vec - 2);
  cout << "vec2 = vec * 4" << endl;
  vec_dump(vec2 = vec * 4);
  cout << "vec2 / 2" << endl;
  vec_dump(vec2 / 2);

  separate();
  cout << "scalar-vec" << endl;
  vec_set(vec);
  cout << "vec: " << endl;
  vec_dump(vec);
  cout << "2 + vec" << endl;
  vec_dump(2 + vec);
  cout << "2 - vec" << endl;
  vec_dump(2 - vec);
  cout << "vec2 = 4 * vec" << endl;
  vec_dump(vec2 = 4 * vec);
  cout << "(vec = 4) / 2" << endl;
  // JADP 17 Maart 2002 vec_dump( (vec = 4) / 2 );
  vec_dump( ( vec.assignToAll(4) ) / 2 );

}

template<typename TF_Vector>
void vec_testcomplex(const TF_Vector&) {

  typedef typename gCppTraits<typename TF_Vector::T_Value>::T_Complex T_Complex;
  typedef gRowVector<T_Complex, DENSE> T_Vector;

  T_Vector vec(10);
  vec_set( vec, T_Complex(0,1), T_Complex(1,1) );
  cout << "vec:" << endl;
  vec_dump(vec);
  cout << "conj(vec):" << endl;
  vec_dump( conj(vec) );
  cout << "real(vec):" << endl;
  vec_dump( real(vec) );
  cout << "imag(vec):" << endl;
  vec_dump( imag(vec) );
  cout << "abs(vec):" << endl;
  vec_dump( abs(vec) );
  cout << "arg(vec)" << endl;
  vec_dump( arg(vec) );
  cout << "norm(vec)" << endl;
  vec_dump( norm(vec) );
  cout << "cos(vec)" << endl;
  vec_dump( cos(vec) );
  cout << "cosh(vec)" << endl;
  vec_dump( cosh(vec) );
  cout << "sin(vec)" << endl;
  vec_dump( sin(vec) );
  cout << "sinh(vec)" << endl;
  vec_dump( sinh(vec) );
  cout << "exp(vec)" << endl;
  vec_dump( exp(vec) );
  cout << "log(vec)" << endl;
  vec_dump( log(vec) );
  cout << "sqrt(vec)" << endl;
  vec_dump( sqrt(vec) );
  cout << "step(real(vec))" << endl;
  vec_dump(step(real(vec)));
  cout << "impulse(real(vec))" << endl;
  vec_dump(impulse(real(vec)));
}

void vec_testdot(void) {

  typedef gRowVector<int, DENSE> T_Rvd;
  typedef gColVector<int, DENSE> T_Cvd;
  typedef gRowVector<int, SLICE> T_Rvs;
  typedef gColVector<int, SLICE> T_Cvs;
  T_Rvd rv(N);
  T_Cvd cv(N);

  cout << "Testing vector-vector dotproduct" << endl;
  separate();
  cout << "DENSE - DENSE" << endl;
  separate();
  vec_set(rv, 1);
  vec_set(cv, 2);
  cout << "rv, cv" << endl;
  vec_dump(rv);
  vec_dump(cv);
  cout << "rv * cv" << endl;
  cout << (rv * cv) << " | " << 5850 << endl;
  cout << "rv * rv.transpose()" << endl;
  cout << (rv * rv.transpose()) << " | " << 5525 << endl;
  cout << "cv.tranpose() * cv" << endl;
  cout << (cv.transpose() * cv) << " | " << 6200 << endl;
  separate();

  separate();
  cout << "SLICE - SLICE" << endl;
  separate();
  T_Rvs rvs = rv.slice(gIndexRange(0, 24, 2));
  T_Cvs cvs = cv.slice(gIndexRange(0, 24, 2));
  cout << "rvs, cvs" << endl;
  vec_dump(rvs);
  vec_dump(cvs);
  cout << "rvs * cvs" << endl;
  cout << (rvs * cvs) << " | " << 3094 << endl;
  cout << "rvs * rvs.transpose()" << endl;
  cout << (rvs * rvs.transpose()) << " | " << 2925 << endl;
  cout << "cvs.tranpose() * cvs" << endl;
  cout << (cvs.transpose() * cvs) << " | " << 3276 << endl;
  separate();

  separate();
  cout << "DENSE - SLICE" << endl;
  separate();
  T_Rvd rvds = rv.subrange(0,12);
  cout << "rvds, cvs" << endl;
  vec_dump(rvds);
  vec_dump(cvs);
  cout << "rvds * cvs" << endl;
  cout << (rvds * cvs) << " | " << 1638 << endl;

  separate();
  cout << "SLICE - DENSE" << endl;
  separate();
  T_Cvd cvds = cv.subrange(0,12);
  cout << "rvs, cvds" << endl;
  vec_dump(rvs);
  vec_dump(cvds);
  cout << "rvs * cvds" << endl;
  cout << (rvs * cvds) << " | " << 1716 << endl;
  separate();

}

template<typename T_Matrix>
void mat_testctor(const T_Matrix&) {

  int r = 1; // keep ref count
  separate();
  cout << "ctors" << endl;
  separate();
  cout << "default ctor (no args)" << endl;
  T_Matrix mat;
  cout << "[mat.rows() x mat.cols()]: " << "[ " << mat.rows() << " x " << mat.cols() <<
    " ] | [ 0 x 0 ]" << endl;

  separate();
  cout << "dtor" << endl;
  mat.~T_Matrix();
  cout << "ctor specifying size (using placement new)" << endl;
  new(&mat) T_Matrix(NR, NC);
  cout << "[mat.rows() x mat.cols()]: " << "[ " << mat.rows() << " x " << mat.cols()
       << " ] | [ " << NR << " x " << NC << " ]" << endl;
  cout << "initialize" << endl;
  mat_set(mat);
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "copy ctor: mat2(mat)" << endl;
  T_Matrix mat2(mat);
  r++;
  cout << "[mat2.rows() x mat2.cols()]: " << "[ " << mat2.rows() << " x " << mat2.cols() <<
    " ] | [ " << NR << " x " << NC << " ]" << endl;
  cout << "mat, mat2" << endl;
  mat_dump(mat);
  mat_dump(mat2);
  cout << "mat2.getRefs(): " << mat2.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "creating and init'ing mat3 (offset 2)" << endl;
  T_Matrix mat3(NR, NC);
  mat_set(mat3, 2);
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "aliasing mat3: mat3.alias(mat):" << endl;
  mat3.alias(mat);
  r++;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat3.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat += 3:" << endl;
  mat += 3;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "deepening mat3: mat3.deepen(): " << endl;
  mat3.deepen();
  r--;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat3.getRefs(): " << mat3.getRefs() << " | " << 1 << endl;
  mat -= 2;
  mat3 += 2;
  cout << "mat -= 2, mat3 += 2" << endl;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);

  separate();
  cout << "assignment: matrix-matrix" << endl;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat3.getRefs(): " << mat3.getRefs() << " | " << 1 << endl;
  cout << "mat3 = mat" << endl;
  mat3 = mat;
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat3.getRefs(): " << mat3.getRefs() << " | " << 1 << endl;
  cout << "assignment: matrix-scalar" << endl;
  cout << "mat3 = 3" << endl;
  // JADP 17 Maart 2002 mat3 = 3;
  mat3.assignToAll(3);
  cout << "mat, mat3" << endl;
  mat_dump(mat);
  mat_dump(mat3);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mat3.getRefs(): " << mat3.getRefs() << " | " << 1 << endl;

  separate();
  cout << "mat: " << endl;
  mat_dump(mat);
  cout << "transpose: mattrans = mat.transpose()" << endl;
  typename T_Matrix::T_Transpose mattrans = mat.transpose();
  r++;
  cout << "[mattrans.rows() x mattrans.cols()]: " << "[ " << mattrans.rows() << " x " << mattrans.cols() <<
    " ] | [ " << NC << " x " << NR << " ]" << endl;
  cout << "mat, mattrans:" << endl;
  mat_dump(mat);
  mat_dump(mattrans);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << "mattrans.getRefs(): " << mattrans.getRefs() << " | " << r << endl;
}

template<typename T_Matrix>
void mat_testsub(const T_Matrix&) {

  int r = 1;
  typedef typename T_Matrix::T_Slice T_Slice;
  typedef typename T_Matrix::T_Rowsub T_Rowsub;
  typedef typename T_Matrix::T_Colsub T_Colsub;

  T_Matrix mat(LNR, LNC);
  mat_set(mat);

  separate();
  cout << "subranges and slices" << endl;
  separate();
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "iterating matrix:" << endl;
  mat_iter_row(mat);
  mat_iter_col(mat);
  cout << "mat.getRefs(): " << mat.getRefs() << " | " << r << endl;
  cout << endl;

  separate();
  cout << "subranging matrix rows: [" << MSRS << ":" << MSRE << ",:]" << endl;
  T_Rowsub matrowsub = mat.extractRows(MSRS, MSRE);
  r++;
  cout << "[matrowsub.rows() x matrowsub.cols()]: " << "[ " << matrowsub.rows() << " x " << matrowsub.cols() <<
    " ] | [ " << MSRE - MSRS + 1 << " x " << mat.cols() << " ]" << endl;
  cout << "matrowsub: " << endl;
  mat_dump(matrowsub);
  cout << "iterating matrowsub: " << endl;
  mat_iter_row(matrowsub);
  mat_iter_col(matrowsub);
  cout << "matrowsub.getRefs(): " << matrowsub.getRefs() << " | " << r << endl;
  cout << endl;

  cout << "subranging matrix colunms: [:," << MSCS << ":" << MSCE << "]" << endl;
  T_Colsub matcolsub = mat.extractCols(MSCS, MSCE);
  r++;
  cout << "[matcolsub.rows() x matcolsub.cols()]: " << "[ " << matcolsub.rows() << " x " << matcolsub.cols() <<
    " ] | [ " << mat.rows() << " x " << MSCE - MSCS + 1 << " ]" << endl;
  cout << "matcolsub: " << endl;
  mat_dump(matcolsub);
  cout << "iterating matcolsub: " << endl;
  mat_iter_row(matcolsub);
  mat_iter_col(matcolsub);
  cout << "matcolsub.getRefs(): " << matcolsub.getRefs() << " | " << r << endl;
  cout << endl;

  cout << "slicing matrix: [" << MSRS << ":" << MSRT << ":" << MSRE << "," <<
    MSCS << ":" << MSCT << ":" << MSCE << "]" << endl;
  T_Slice matslice = mat.slice(gIndexRange(MSRS, MSRE, MSRT),
             gIndexRange(MSCS, MSCE, MSCT));
  r++;
  cout << "[matslice.rows() x matslice.cols()]: " << "[ " << matslice.rows() << " x " << matslice.cols() <<
    " ] | [ " <<
    ceil(static_cast<double>(MSRE - MSRS + 1) / static_cast<double>(MSRT)) << " x " <<
    ceil(static_cast<double>(MSCE - MSCS + 1) / static_cast<double>(MSCT)) << " ]" << endl;
  cout << "matslice:" << endl;
  mat_dump(matslice);
  cout << "iterating matslice" << endl;
  mat_iter_row(matslice);
  mat_iter_col(matslice);
  cout << "matslice.getRefs(): " << matslice.getRefs() << " | " << r << endl;
  cout << endl;
}

template<typename T_Matrix>
void mat_testopers(const T_Matrix&) {

  typedef typename T_Matrix::T_Value T_Value;
  typedef gRowVector<T_Value, DENSE> T_Vector;

  T_Matrix mat(NR, NC);
  T_Matrix mat2(NR, NC);

  separate();
  cout << "Matrix operators and functions" << endl;
  separate();
  cout << "Unary minus" << endl;
  mat_set(mat);
  cout << "mat, -mat" << endl;
  mat_dump(mat);
  mat_dump(-mat);
  cout << "Ones" << endl;
  mat_dump(ones<T_Value>(NR, NC));
  cout << "Zeros" << endl;
  mat_dump(zeros<T_Value>(NR, NC));
  cout << "Identity matrix" << endl;
  mat_dump(eye<T_Value>(NR, NC));
  mat_dump(eye<T_Value>(NC, NR));
  mat_dump(eye<T_Value>(NC));
  cout << "Diagonal matrix from vector" << endl;
  T_Vector vec(NC);
  vec_set(vec,1);
  mat_dump(diag(vec));
  mat_dump(diag(vec.transpose()));

  T_Vector vec2(NR);
  vec_set(vec);
  vec_set(vec2);

  cout << "Outer product" << endl;
  cout << "vec, vec2" << endl;
  vec_dump(vec);
  vec_dump(vec2);
  cout << "vec2.transpose() * vec" << endl;
  mat_dump(vec2.transpose() * vec);
  cout << "Equality and inequality operator" << endl;
  mat_set(mat);
  mat_set(mat2);
  cout << "mat, mat2" << endl;
  mat_dump(mat);
  mat_dump(mat2);
  cout << (mat == mat2) << " " << (mat != mat2) << " | " << "1 0" << endl;
  mat2 += 1;
  cout << "mat, mat2" << endl;
  mat_dump(mat);
  mat_dump(mat2);
  cout << (mat == mat2) << " " << (mat != mat2) << " | " << "0 1" << endl;
  cout << "mat + mat2" << endl;
  mat_dump(mat + mat2);
  cout << "mat - mat2" << endl;
  mat_dump(mat - mat2);
  cout << "mat * mat2.transpose()" << endl;
  mat_dump(mat * mat2.transpose());
  cout << "(mat * mat2.transpose()).transpose() == (mat2 * mat.transpose())" << endl;
  cout << ( (mat * mat2.transpose()).transpose() == (mat2 * mat.transpose()) ) << endl;
  cout << "mat += mat2" << endl;
  mat += mat2;
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "mat -= mat2" << endl;
  mat -= mat2;
  mat_dump(mat);
  cout << "matrix - colvector multiply" << endl;
  mat_dump(mat);
  vec.resize(mat.cols());
  vec_set(vec, 1);
  vec_dump(vec);
  vec_dump(mat * vec.transpose());
  cout << "rowvector - matrix multiply" << endl;
  mat_dump(mat);
  vec.resize(mat.rows());
  vec_set(vec, 1);
  vec_dump(vec);
  vec_dump(vec * mat);
  mat.resize(2,2);
  mat_set(mat, 1);
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "mat + 2" << endl;
  mat_dump(mat + 2);
  cout << "mat - 2" << endl;
  mat_dump(mat - 2);
  cout << "mat * 2" << endl;
  mat_dump(mat * 2);
  cout << "mat / 2" << endl;
  mat_dump(mat / 2);
  cout << "2 + mat" << endl;
  mat_dump(2 + mat);
  cout << "2 - mat" << endl;
  mat_dump(2 - mat);
  cout << "2 * mat" << endl;
  mat_dump(2 * mat);
  cout << "10 / mat" << endl;
  mat_dump(10 / mat);
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "mat += 8" << endl;
  mat += 8;
  mat_dump(mat);
  cout << "mat -= 6" << endl;
  mat -= 6;
  mat_dump(mat);
  cout << "mat *= 4" << endl;
  mat *= 4;
  mat_dump(mat);
  cout << "mat /= 3" << endl;
  mat /= 3;
  mat_dump(mat);
}

template<typename TF_Matrix>
void mat_testcomplex(const TF_Matrix&) {

  typedef typename gCppTraits<typename TF_Matrix::T_Value>::T_Complex T_Complex;
  typedef gMatrix<T_Complex, typename TF_Matrix::T_Storage> T_Matrix;

  T_Matrix mat(3, 4);
  mat_set(mat, T_Complex(0,1), T_Complex(1,1));
  cout << "mat:" << endl;
  mat_dump(mat);
  cout << "conj(mat):" << endl;
  mat_dump(conj(mat));
  cout << "real(mat):" << endl;
  mat_dump(real(mat));
  cout << "imag(mat):" << endl;
  mat_dump(imag(mat));
  cout << "abs(mat):" << endl;
  mat_dump(abs(mat));
  cout << "arg(mat)" << endl;
  mat_dump(arg(mat));
  cout << "norm(mat)" << endl;
  mat_dump(norm(mat));
  cout << "cos(mat)" << endl;
  mat_dump(cos(mat));
  cout << "cosh(mat)" << endl;
  mat_dump(cosh(mat));
  cout << "sin(mat)" << endl;
  mat_dump(sin(mat));
  cout << "sinh(mat)" << endl;
  mat_dump(sinh(mat));
  cout << "exp(mat)" << endl;
  mat_dump(exp(mat));
  cout << "log(mat)" << endl;
  mat_dump(log(mat));
  cout << "sqrt(mat)" << endl;
  mat_dump(sqrt(mat));
  //----------------------lude not yet in out.reference  22/11/2001
  cout << "step(real(mat))" << endl;
  mat_dump(step(real(mat)));
  cout << "impulse(real(mat))" << endl;
  mat_dump(impulse(real(mat)));
  //----------------------lude not yet in out.reference  22/11/2001
}

void testlapack(void) {

  typedef double T_V;
  typedef gMatrix<T_V, ROW_DENSE> T_Matrix;

  T_Matrix mat(3,3);
  mat(0,0) = 2; mat(0,1) = 1; mat(0,2) = 1;
  mat(1,0) = 1; mat(1,1) = 2; mat(1,2) = 1;
  mat(2,0) = 1; mat(2,1) = 1; mat(2,2) = 2;

  T_Matrix a;
  a = mat;
  cout << "a:" << endl;
  mat_dump(a);
  cout << "inv(a):" << endl;
  int fail;
  mat_dump(inv(a, fail));
  gEigenData<complex<T_V>, complex<T_V> > eigs = eig(a);
  cout << "eig(a).vectors())" << endl;
  mat_dump(eigs.vectors());
  cout << "eig(a).values())" << endl;
  mat_dump(eigs.values());
}

void testLU(void) {

  cout << endl << "======================================================" << endl;
  cout << "Testing Lapack LU" << endl;

  int fail;

  gMatrix<float, COL_DENSE> A(2, 2);
  A(0, 0) = 4;
  A(0, 1) = 3;
  A(1, 0) = 6;
  A(1, 1) = 3;
  cout << "A = " << endl;
  mat_dump(A);

  gMatrix<float> LU(2, 2);
//   lu(A, LU);

//   cout << "LU = " << endl;
//   mat_dump(LU);

  A.resize(3, 3);
  A(0, 0) = 3;
  A(0, 1) = 6;
  A(0, 2) = -9;
  A(1, 0) = 2;
  A(1, 1) = 5;
  A(1, 2) = -3;
  A(2, 0) = -4;
  A(2, 1) = 1;
  A(2, 2) = 10;

  cout << "A = " << endl;
  mat_dump(A);

  LU.resize(3, 3);
  gRowVector<int> P(3);

  lu_factor(A, LU, P, fail);

  cout << "LU = " << endl;
  mat_dump(LU);

  cout << "P = ";
  vec_dump(P);

  gMatrix<float> L(3, 3);
  gMatrix<float> U(3, 3);
  L.assignToAll(0);
  U.assignToAll(0);

  U(0, 0) = LU(0, 0);
  U(0, 1) = LU(0, 1);
  U(0, 2) = LU(0, 2);
  U(1, 1) = LU(1, 1);
  U(1, 2) = LU(1, 2);
  U(2, 2) = LU(2, 2);

  L(0, 0) = 1;
  L(1, 1) = 1;
  L(2, 2) = 1;

  L(1, 0) = LU(1, 0);
  L(2, 0) = LU(2, 0);
  L(2, 1) = LU(2, 1);

  cout << "L = " << endl;
  mat_dump(L);
  cout << "U = " << endl;
  mat_dump(U);

  gMatrix<float> B( L*U );
  cout << "L*U = " << endl;
  mat_dump(B);

  gMatrix<float,COL_DENSE> PA(A);
  PA.deepen();
  lu_permute(PA, P);
  cout << "P*A == L*U " << endl;
  mat_dump(PA);

  gRowVector<int> revP(3);
  for (int i = 0; i < 3; i++) {
    revP[ P[i] ] = i;
  }

  B = L*U;
  lu_permute(B, revP);
  cout << "P^-1*L*U == A:" << endl;
  mat_dump(B);


  cout << "======================================================" << endl;
  cout << "Testing Lapack LU Solve: " << endl;
  gColVector<float> answ1(3);
  gColVector<float> answ2(3);
  answ1.set(1,1);
  answ2.set(10,-1);
  gColVector<float> b1 = A*answ1;
  gColVector<float> b2 = A*answ2;

  cout << "A = " << endl;  mat_dump(A);
  cout << "rhs1 = ";  vec_dump(b1);
  cout << "rhs2 = ";  vec_dump(b2);
  cout << "true x1 = ";  vec_dump(answ1);
  cout << "true x2 = ";  vec_dump(answ2);
  cout << endl;
  gColVector<float> x(3);
  x.assignToAll(0);

  gMatrix<float> bb(3,2);
  gMatrix<float> xx(3,2);

  bb(0,0) = b1[0];
  bb(1,0) = b1[1];
  bb(2,0) = b1[2];
  bb(0,1) = b2[0];
  bb(1,1) = b2[1];
  bb(2,1) = b2[2];

  lu_solve(A, b1, x, fail);
  cout << "Solution 1 = "; vec_dump(x);
  lu_solve(A, b2, x, fail);
  cout << "Solution 2 = "; vec_dump(x);

  lu_solve(LU, P, b1, x, fail);
  cout << "Solution with precalculated LU 1: "; vec_dump(x);
  lu_solve(LU, P, b2, x, fail);
  cout << "Solution with precalculated LU 2: "; vec_dump(x);

  cout << "Solution via multiple rhs (in columns): " << endl;
  lu_solve(A, bb, xx, fail);
  mat_dump(xx);

  cout << "Solution via multiple rhs (in columns) with precalculated LU factorization: " << endl;
  xx.assignToAll(0.0);
  lu_solve(LU, P, bb, xx, fail);
  mat_dump(xx);

  cout << "Solution via multiple rhs (in columns) with precalculated LU factorization: " << endl;
  cout << "Solution with permuted rhs, identity permutation P = [0,1,2]: " << endl;
  cout << " => solve PAx = LUx = Pb" << endl;
  xx.assignToAll(0.0);
  // permute rhs
  gMatrix<float> Pbb(bb);
  Pbb.deepen();
  lu_permute(Pbb, P);
  cout << "P*b = " << endl;
  mat_dump(Pbb);
  // Identiy permutation eyeP = I
  gRowVector<int>  eyeP(3);
  eyeP.set(0,1);

  lu_solve(LU, eyeP, Pbb, xx, fail);
  cout << "Solution = " << endl;
  mat_dump(xx);

  cout << endl << "======================================================" << endl;
  cout << "Testing Lapack cholesky" << endl;

  L.resize(3,3); L.assignToAll(0.0);
  L(0,0) = 1;
  L(1,0) = 2; L(1,1) = 4;
  L(2,0) = 3; L(2,1) = 5; L(2,2) = 6;

  A = L*L.transpose();
  gColVector<float> ans(3);
  ans.set(1,1);
  gColVector<float> b(A*ans);

  cout << "Init A: " << A;
  cout << "Init L: " << L;
  cout << "Init x: " << x;
  cout << "Init b: " << b;

  cout << endl << "Cholesky factor test: " << endl;
  gMatrix<float> chol( A.rows(), A.cols() );
  cholesky_factor(A, chol, fail);
  cout << "L: " << chol;

  cout << endl << "Cholesky solve single rhs test: " << endl;
  x.resize(3);
  x.assignToAll(0.0);
  cholesky_solve(A, b, x, fail);
  cout << "x: " << x;

  cout << endl << "Cholesky solve multiple rhs test: " << endl;
  B.resize(3,2);
  B.col(0) = B.col(1) = b;
  gMatrix<float> X;
  cholesky_solve(A, B, X, fail);
  cout << "X: " << X;

  cout << endl << "Cholesky solve with precomputed L: " << endl;
  fortran::fortran_array<fortran::real> fL(9); // keep in Fortran structure
  cholesky_factor(A, fL, fail);

  x.assignToAll(0.0);
  cholesky_solve(fL, b, x, fail);
  cout << "x: " << x;

  cout << endl << "Cholesky solve multiple rhs with precomputed L: " << endl;
  X.assignToAll(0.0);
  cholesky_solve(fL, B, X, fail);
  cout << "X: " << X;

  cout << "Done, Cholesky tests." << endl;
} // testLU




int main(void) {

#ifdef TEST_SPECIAL


#else // TEST_SPECIAL

  typedef const gRowVector<int, DENSE> T_Vec;
  typedef gRowVector<int, DENSE>::T_Slice T_Slice;
  typedef T_Slice::T_Iterator T_Iter;

  separate();
  cout << "Testing RowVector" << endl;
  vec_testctor(gRowVector<int, DENSE>(gLinearNoInit()));
  vec_testsub(gRowVector<int, DENSE>(gLinearNoInit()));
  vec_testmemopers(gRowVector<int, DENSE>(gLinearNoInit()));
  vec_testopers(gRowVector<int, DENSE>(gLinearNoInit()));
  //----------------------lude not yet in out.reference  22/11/2001
  vec_testcomplex(gRowVector<float, DENSE>(gLinearNoInit()));
  //----------------------lude not yet in out.reference  22/11/2001

  separate();
  cout << "Testing ColVector" << endl;
  vec_testctor(gColVector<int, DENSE>(gLinearNoInit()));
  vec_testsub(gColVector<int, DENSE>(gLinearNoInit()));
  vec_testmemopers(gColVector<int, DENSE>(gLinearNoInit()));
  vec_testopers(gColVector<int, DENSE>(gLinearNoInit()));

  cout << "Testing RowVector and ColVector" << endl;
  vec_testdot();

  separate();
  cout << "Testing Matrix<ROW_DENSE>" << endl;
  mat_testctor(gMatrix<int, ROW_DENSE>(gLinearNoInit()));
  mat_testsub(gMatrix<int, ROW_DENSE>(gLinearNoInit()));
  mat_testopers(gMatrix<int, ROW_DENSE>(gLinearNoInit()));
  mat_testcomplex(gMatrix<float, ROW_DENSE>(gLinearNoInit()));

  cout << "Testing Matrix<COL_DENSE>" << endl;
  mat_testctor(gMatrix<int, COL_DENSE>(gLinearNoInit()));
  mat_testsub(gMatrix<int, COL_DENSE>(gLinearNoInit()));
  mat_testopers(gMatrix<int, COL_DENSE>(gLinearNoInit()));
  mat_testcomplex(gMatrix<float, COL_DENSE>(gLinearNoInit()));

  cout << "Testing Lapack" << endl;
  testlapack();

  testLU();
  cout << "Done." << endl;

#endif
}
