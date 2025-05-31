/*
  $Id$
*/

 /*****************************************************************************
 *               |                                                            *
 * Author        | Ludwig Schwardt                                            *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 2001/11/08                                                 *
 * Description   | Optimized outer products and quadratic forms               *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

#ifdef __ARM_NEON__
#include <arm_neon.h>
#endif

// patrec headers
#include "prlite_gendefs.hpp"  // REAL
#include "prlite_genvec.hpp"  // Vector
#include "prlite_error.hpp"  // FAIL
#include "prlite_genmat.hpp"  // Matrix

namespace prlite{

// Copies lower triangle to upper triangle, thereby completing symmetric matrix.
template<typename Container>
inline void copyLowerToUpper(Container&) {
} // copyLowerToUpper

// Extend lower triangular matrix to full symmetric one
template<typename T>
void copyLowerToUpper(RowMatrix<T>& A) {
  int sz( A.rows() );
#ifdef CHK_IDX
  if (A.cols() != sz) {
    PRLITE_FAIL("copyLowerToUpper(Matrix) : matrix not square "
         << "(" << A.rows() << "x" << A.cols() << " instead)");
  } // if
#endif // CHK_IDX
  typename RowMatrix<T>::T_RowIterator outrowtoitr( A.rowBegin() );
  typename RowMatrix<T>::T_Row::T_Iterator outcoltoitr;
  typename RowMatrix<T>::T_ConstColIterator outcolfromitr( A.colBegin() );
  typename RowMatrix<T>::T_Col::T_ConstIterator outrowfromitr;
  int cols;
  int rows(sz);
  while (--rows) {
    cols = rows;
    outcoltoitr = (*outrowtoitr++).begin() + (sz - cols);
    outrowfromitr = (*outcolfromitr++).begin() + (sz - cols);
    while (cols--) {
      *outcoltoitr++ = *outrowfromitr++;
    } // while
  } // while
} // copyLowerToUpper

// Extend lower triangular matrix to full symmetric one
template<typename T, typename T_STORAGE>
void copyLowerToUpper(gLinear::gMatrix<T, T_STORAGE>& A) {
  int sz( A.rows() );
#ifdef CHK_IDX
  if (A.cols() != sz) {
    PRLITE_FAIL("copyLowerToUpper(Matrix) : matrix not square "
         << "(" << A.rows() << "x" << A.cols() << " instead)");
  } // if
#endif // CHK_IDX
  typename gLinear::gMatrix<T,T_STORAGE>::T_RowIterator outrowtoitr( A.rowBegin() );
  typename gLinear::gMatrix<T,T_STORAGE>::T_Row::T_Iterator outcoltoitr;
  typename gLinear::gMatrix<T,T_STORAGE>::T_ConstColIterator outcolfromitr( A.colBegin() );
  typename gLinear::gMatrix<T,T_STORAGE>::T_Col::T_ConstIterator outrowfromitr;
  int cols;
  int rows(sz);
  while (--rows) {
    cols = rows;
    outcoltoitr = (*outrowtoitr++).begin() + (sz - cols);
    outrowfromitr = (*outcolfromitr++).begin() + (sz - cols);
    while (cols--) {
      *outcoltoitr++ = *outrowfromitr++;
    } // while
  } // while
} // copyLowerToUpper

} // namespace prlite
