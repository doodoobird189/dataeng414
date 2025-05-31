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

#ifndef PRLITE_OUTERQUAD_L2U_HPP
#define PRLITE_OUTERQUAD_L2U_HPP

// patrec headers
#include "prlite_gendefs.hpp"  // REAL
#include "prlite_genvec.hpp"  // Vector
#include "prlite_genmat.hpp"  // Matrix

namespace prlite{

/**
 * Copies lower triangle to upper triangle, thereby completing symmetric matrix.
 * This extends a lower triangular matrix to full symmetric form. Useful to process
 * the output of sumOuterProdSym. Does nothing by default (i.e. if input not matrix).
 * @param A Matrix to be extended (should be lower triangular)
 */
template<typename Container>
void copyLowerToUpper(Container& A);

// Extend lower triangular matrix to full symmetric one
template<typename T>
void copyLowerToUpper(RowMatrix<T>& A);

// Extend lower triangular matrix to full symmetric one
template<typename T, typename T_STORAGE>
void copyLowerToUpper(gLinear::gMatrix<T,T_STORAGE>& A);

} // namespace prlite

#include "prlite_outerquad_l2u.tcc"

#endif // PRLITE_OUTERQUAD_L2U_HPP
