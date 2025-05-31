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
 * $Id$
 */

// standard headers
#include <cstddef>  // NULL
#include <new>  // placement new
#include <iostream>  // cout, endl


GLINEAR_NAMESPACE_BEGIN

// no init ctor.
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(gLinearNoInit) :
  stor(NULL) {
}

// dflt ctor. [0 x 0]. allocate dflt storage. incr ref cnt.
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(void) :
  nRows(0),
  nCols(0) {
  stor = new T_StorageElem(0);
  stor->addref();
  startAddr = stor->getStartAddr();
}

// ctor specifying size. [rows x cols]. allocate storage size (rows*cols).
// inc ref cnt.
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(T_Index _rows,
                                          T_Index _cols) :
  nRows(_rows),
  nCols(_cols) {
  stor = new T_StorageElem(_rows * _cols);
  stor->addref();
  startAddr = stor->getStartAddr();
}

// copy ctor. shallow copy. clone arg Matrix. inc ref cnt.
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(const gMatrixMemBase<TF_Matrix>& mat) :
  nRows(mat.nRows),
  nCols(mat.nCols),
  stor(mat.stor),
  startAddr(mat.startAddr) {
  stor->addref();
}

// ctor with supplied values (init ctor)
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(T_Index new_nRows,
                                          T_Index new_nCols,
                                          const T_StorageElem* new_stor,
                                          T_ConstPointer new_startAddr) :
  nRows(new_nRows),
  nCols(new_nCols),
  stor( const_cast<T_StorageElem*>(new_stor) ),
  startAddr( const_cast<T_Pointer>(new_startAddr) ) {
  stor->addref();
}

template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(T_ConstPointer addr,
                                          T_Index rows, T_Index cols) :
  nRows(rows),
  nCols(cols),
  stor(new T_StorageElem(addr, rows*cols)),
  startAddr(const_cast<T_Pointer>(stor->getStartAddr())) {
  stor->addref();
}

template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::gMatrixMemBase(T_ConstPointer addr,
                                          T_Index rows, T_Index cols,
                                          const matrix_no_copy_tag&) :
  nRows(rows),
  nCols(cols),
  stor(new T_StorageElem()),
  startAddr(const_cast<T_Pointer>(addr)) {
  stor->setExisting(addr, rows*cols);
}

// dtor. decr ref cnt. delete storage.
template<typename TF_Matrix>
gMatrixMemBase<TF_Matrix>::~gMatrixMemBase(void) {
  if (stor) {
    stor->rmref();
  } // if
}

// configure with supplied values using init ctor.
template<typename TF_Matrix>
void gMatrixMemBase<TF_Matrix>::configure(T_Index new_nRows,
                                          T_Index new_nCols,
                                          const T_StorageElem* new_stor,
                                          T_ConstPointer new_startAddr) {
  // destruct and use placement-new to call init ctor to configure.
  this->~gMatrixMemBase();
  new(this) gMatrixMemBase(new_nRows, new_nCols, new_stor, new_startAddr);
}

// dump info
template<typename TF_Matrix>
void gMatrixMemBase<TF_Matrix>::specialInfo(int i) const {
  std::cout << "gMatrixMemBase<T_Matrix>" << std::endl;
  std::cout << "  [nRows x nCols]: " << "[ " << nRows << " x " << nCols << " ]"
       << " stor: " << stor
       << " startAddr: " << startAddr
       << std::endl;
  std::cout << " gMatrixMemBase<T_Matrix>::";
  stor->info(i);
}

// return number of rows
template<typename TF_Matrix>
typename gMatrixMemBase<TF_Matrix>::T_Index
gMatrixMemBase<TF_Matrix>::specialRows(void) const {
  return nRows;
}

// return number of columns
template<typename TF_Matrix>
typename gMatrixMemBase<TF_Matrix>::T_Index
gMatrixMemBase<TF_Matrix>::specialCols(void) const {
  return nCols;
}

// resize
template<typename TF_Matrix>
bool gMatrixMemBase<TF_Matrix>::specialResize(T_Index new_nRows,
                                              T_Index new_nCols) {
  nRows = new_nRows;
  nCols = new_nCols;
  stor->setSize(new_nRows * new_nCols);
  startAddr = stor->getStartAddr();
  return true;
}

// return number of refereces to the storage object
template<typename TF_Matrix>
typename gMatrixMemBase<TF_Matrix>::T_Index
gMatrixMemBase<TF_Matrix>::specialGetRefs(void) const {
  return stor->getRefs();
}

// ensure no aliases exist
template<typename TF_Matrix>
inline void
gMatrixMemBase<TF_Matrix>::specialDeepen(void) const {
  // save starting index.
  T_Index delta = startAddr - stor->getStartAddr();
  // duplicate stor
  T_StorageElem* tmp = stor->duplicate();
  // decrement ref on old stor, reassign stor and inc ref on new stor
  stor->rmref();
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->stor = tmp;
  stor->addref();
  // restore starting index
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->startAddr = stor->getStartAddr() + delta;
}

// const alias
template<typename TF_Matrix>
inline void
gMatrixMemBase<TF_Matrix>::specialAlias(const gMatrixMemBase<TF_Matrix>& vec) const {
  // set size
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->nRows = vec.nRows;
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->nCols = vec.nCols;
  // decr ref on old stor, reassign stor and incr ref on new stor.
  // Note: this is not atomic.
  stor->rmref();
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->stor = vec.stor;
  stor->addref();
  // set start address
  GLINEAR_CONST_OVERRIDE_THIS(gMatrixMemBase<TF_Matrix>)->startAddr = vec.startAddr;
}

// alias
template<typename TF_Matrix>
inline void
gMatrixMemBase<TF_Matrix>::specialAlias(gMatrixMemBase<TF_Matrix>& vec) {
  // set size
  nRows = vec.nRows;
  nCols = vec.nCols;
  // decr ref on old stor, reassign stor and incr ref on new stor.
  // Note: this is not atomic.
  stor->rmref();
  stor = vec.stor;
  stor->addref();
  // set start address
  startAddr = vec.startAddr;
}

// return const stor
template<typename TF_Matrix>
inline const typename gMatrixMemBase<TF_Matrix>::T_StorageElem*
gMatrixMemBase<TF_Matrix>::getStor(void) const {
  return stor;
}

// return stor
template<typename TF_Matrix>
inline typename gMatrixMemBase<TF_Matrix>::T_StorageElem*
gMatrixMemBase<TF_Matrix>::getStor(void) {
  return stor;
}

// return const startAddr
template<typename TF_Matrix>
inline typename gMatrixMemBase<TF_Matrix>::T_ConstPointer
gMatrixMemBase<TF_Matrix>::getStartAddr(void) const {
  return startAddr;
}

// return startAddr
template<typename TF_Matrix>
inline typename gMatrixMemBase<TF_Matrix>::T_Pointer
gMatrixMemBase<TF_Matrix>::getStartAddr(void) {
  return startAddr;
}

GLINEAR_NAMESPACE_END
