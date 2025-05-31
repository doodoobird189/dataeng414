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
#include <algorithm>  // transform
#include <iostream>  // ostream, ios, endl, istream


GLINEAR_NAMESPACE_BEGIN

// rowvector unary minus
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
operator-(const gRowVector<T, T_Stor>& vec) {
  return o__vecUnaryMinus(vec);
}

// colvector unary minus
template<typename T, typename T_Stor>
typename gColVector<T, T_Stor>::T_Dense
operator-(const gColVector<T, T_Stor>& vec) {
  return o__vecUnaryMinus(vec);
}

// generic unary minus -- support
template<typename T_Vec>
typename T_Vec::T_Dense
o__vecUnaryMinus(const T_Vec& vec) {
  typedef typename T_Vec::T_Dense T_Result;
  T_Result result(vec);
  result.deepen();
  typename T_Result::T_ConstIterator iterend = result.end();
  for (typename T_Result::T_Iterator iter = result.begin(); iter != iterend; iter++) {
    *iter = - (*iter);
  }
  return result;
}

// vector equality
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator==(const gRowVector<T1, lhstortype>& vecl,
           const gRowVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator==(const gRowVector<T, lhstortype>&, const gRowVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    if ( vecl[indx] != vecr[indx] ) {
      return false;
    }
  }
  return true;
}

// vector inequality
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator!=(const gRowVector<T1, lhstortype>& vecl,
           const gRowVector<T2, rhstortype>& vecr) {
  return !(vecl == vecr);
}

// vector less than
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator<(const gRowVector<T1, lhstortype>& vecl,
          const gRowVector<T2, rhstortype>& vecr) {

  if ( vecl.size() < vecr.size() ) { return true; } // if
  if ( vecl.size() > vecr.size() ) { return false; } // if

  typedef gRowVector<T1, lhstortype> T_Vleft;
  typedef gRowVector<T2, rhstortype> T_Vright;
  const typename T_Vleft::T_ConstIterator veclend = vecl.end();
  typename T_Vleft::T_ConstIterator vecliter = vecl.begin();
  typename T_Vright::T_ConstIterator vecriter = vecr.begin();

  while (vecliter != veclend) {
    if (*vecliter < *vecriter) { return true; } // if
    if (*vecliter++ > *vecriter++) { return false; } // if
  } // while

  return false;
} // operator<

// gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> addition
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gRowVector<T1, lhstortype>& vecl,
          const gRowVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator+(const gRowVector<T, lhstortype>&, const gRowVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> tmp;
  tmp = vecl;
  tmp +=vecr;
  return tmp;
}

// gRowVector<T1, lhstortype> - gRowVector<T2, rhstortype> subtraction
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gRowVector<T1, lhstortype>& vecl,
          const gRowVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator-(const gRowVector<T, lhstortype>&, const gRowVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> tmp;
  tmp = vecl;
  tmp -= vecr;
  return tmp;
}

// vector equality
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator==(const gColVector<T1, lhstortype>& vecl,
           const gColVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator==(const gColVector<T, lhstortype>&, const gColVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    if ( vecl[indx] != vecr[indx] ) {
      return false;
    }
  }
  return true;
}

// vector inequality
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator!=(const gColVector<T1, lhstortype>& vecl,
           const gColVector<T2, rhstortype>& vecr) {
  return !(vecl == vecr);
}

// vector lessthan
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
bool
operator<(const gColVector<T1, lhstortype>& vecl,
          const gColVector<T2, rhstortype>& vecr) {

  if ( vecl.size() < vecr.size() ) { return true; } // if
  if ( vecl.size() > vecr.size() ) { return false; } // if

  typedef gColVector<T1, lhstortype> T_Vleft;
  typedef gColVector<T2, rhstortype> T_Vright;
  const typename T_Vleft::T_ConstIterator veclend = vecl.end();
  typename T_Vleft::T_ConstIterator vecliter = vecl.begin();
  typename T_Vright::T_ConstIterator vecriter = vecr.begin();

  while (vecliter != veclend) {
    if (*vecliter < *vecriter) { return true; } // if
    if (*vecliter++ > *vecriter++) { return false; } // if
  } // while

  return false;
} // operator<

// gColVector<T1, stortype> - gColVector<T2, stortype> addition
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gColVector<T1, lhstortype>& vecl,
          const gColVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator+(const gColVector<T, lhstortype>&, const gColVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> tmp;
  tmp = vecl;
  tmp +=vecr;
  return tmp;
}

// gColVector<T1, stortype> - gColVector<T2, stortype> subtraction
template<typename T1, typename T2, typename lhstortype, typename rhstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gColVector<T1, lhstortype>& vecl,
          const gColVector<T2, rhstortype>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator-(const gColVector<T, lhstortype>&, const gColVector<T, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> tmp;
  tmp = vecl;
  tmp -= vecr;
  return tmp;
}

// gRowVector<T1, DENSE> - gColVector<T2, SLICE> multiply (dotpoduct)
template<typename T1, typename T2, typename T_LHStor, typename T_RHStor>
typename gNumPrec2<T1, T2>::T_Highest operator*(const gRowVector<T1, T_LHStor>& vecl,
                                                const gColVector<T2, T_RHStor>& vecr) {
#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator*(const gRowVector<T, DENSE>, const gColVector<T, DENSE>&)",
    "Size mismatch",
    vecl.size(), vecr.size() );
  }
#endif // CHK_IDX
  typedef gRowVector<T1, T_LHStor> T_Vleft;
  typedef gColVector<T2, T_RHStor> T_Vright;
  const typename T_Vleft::T_ConstIterator veclend = vecl.end();
  typename T_Vleft::T_ConstIterator vecliter = vecl.begin();
  typename T_Vright::T_ConstIterator vecriter = vecr.begin();
  typename gNumPrec2<T1, T2>::T_Highest result = 0;
  while (vecliter != veclend) {
    result += (*vecliter++) * (*vecriter++);
  }
  return result;
}

//----------------------jadp start of mods  13/7/2001

// gRowVector<T1, T_Stor> == T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator==(const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] == valr);
  } // while
  return vecret;
} // operator==

// gRowVector<T1, T_Stor> != T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator!=(const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] != valr);
  } // while
  return vecret;
} // operator!=

// gRowVector<T1, T_Stor> > T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator> (const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] > valr);
  } // while
  return vecret;
} // operator>

// gRowVector<T1, T_Stor> < T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator< (const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] < valr);
  } // while
  return vecret;
} // operator<

// gRowVector<T1, T_Stor> >= T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator>=(const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] >= valr);
  } // while
  return vecret;
} // operator>=

// gRowVector<T1, T_Stor> <= T2
template<typename T1, typename T2, typename T_Stor>
gRowVector<bool, DENSE>
operator<=(const gRowVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] <= valr);
  } // while
  return vecret;
} // operator<=

//----------------------jadp end of mods  13/7/2001

// gRowVector<T1, stortype> + T2
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gRowVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return (result += scalar);
}

// gRowVector<T1, stortype> - T2
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gRowVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result -= scalar;
}

// gRowVector<T1, stortype> * T2
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const gRowVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result *= scalar;
}

// gRowVector<T1, stortype> / T2
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const gRowVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result /= scalar;
}


// T1 + gRowVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const T1& scalar,
          const gRowVector<T2, vecstortype>& vec) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result += scalar;
}

// T1 - gRowVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const T1& scalar,
          const gRowVector<T2, vecstortype>& vec) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  typedef typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_ConstIterator T_ConstIterator;
  typedef typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_Iterator T_Iterator;
  const T_ConstIterator dstend = result.end();
  for (T_Iterator dst = result.begin(); dst != dstend; dst++) {
    *dst = scalar - *dst;
  }
  return result;
}

// T1 * gRowVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const T1& scalar,
          const gRowVector<T2, vecstortype>& vec) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result *= scalar;
}

// T1 / gRowVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const T1& scalar,
          const gRowVector<T2, vecstortype>& vec) {
  gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  typedef typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_ConstIterator T_ConstIterator;
  typedef typename gRowVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_Iterator T_Iterator;
  const T_ConstIterator dstend = result.end();
  for (T_Iterator dst = result.begin(); dst != dstend; dst++) {
    *dst = scalar / *dst;
  }
  return result;
}

//----------------------jadp start of mods  13/7/2001

// gColVector<T1, T_Stor> == T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator==(const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] == valr);
  } // while
  return vecret;
} // operator==

// gColVector<T1, T_Stor> != T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator!=(const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] != valr);
  } // while
  return vecret;
} // operator!=

// gColVector<T1, T_Stor> > T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator> (const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] > valr);
  } // while
  return vecret;
} // operator>

// gColVector<T1, T_Stor> < T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator< (const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] < valr);
  } // while
  return vecret;
} // operator<

// gColVector<T1, T_Stor> >= T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator>=(const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] >= valr);
  } // while
  return vecret;
} // operator>=

// gColVector<T1, T_Stor> <= T2
template<typename T1, typename T2, typename T_Stor>
gColVector<bool, DENSE>
operator<=(const gColVector<T1, T_Stor>& vecl,
           const T2& valr) {
  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] <= valr);
  } // while
  return vecret;
} // operator<=

//----------------------jadp end of mods  13/7/2001

// gColVector<T1, stortype> + T2
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const gColVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result += scalar;
}

// gColVector<T1, stortype> - T2
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const gColVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result -= scalar;
}

// gColVector<T1, stortype> * T2
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const gColVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result *= scalar;
}

// gColVector<T1, stortype> / T2
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const gColVector<T1, vecstortype>& vec,
          const T2& scalar) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result /= scalar;
}

// T1 + gColVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator+(const T1& scalar,
          const gColVector<T2, vecstortype>& vec) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result += scalar;
}

// T1 - gColVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator-(const T1& scalar,
          const gColVector<T2, vecstortype>& vec) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  typedef typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_ConstIterator T_ConstIterator;
  typedef typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_Iterator T_Iterator;
  const T_ConstIterator dstend = result.end();
  for (T_Iterator dst = result.begin(); dst != dstend; dst++) {
    *dst = scalar - *dst;
  }
  return result;
}

// T1 * gColVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator*(const T1& scalar,
          const gColVector<T2, vecstortype>& vec) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  return result *= scalar;
}

// T1 / gColVector<T2, stortype>
template<typename T1, typename T2, typename vecstortype>
gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>
operator/(const T1& scalar,
          const gColVector<T2, vecstortype>& vec) {
  gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE> result( vec.size() );
  result = vec;
  typedef typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_ConstIterator T_ConstIterator;
  typedef typename gColVector<typename gNumPrec2<T1, T2>::T_Highest, DENSE>::T_Iterator T_Iterator;
  const T_ConstIterator dstend = result.end();
  for (T_Iterator dst = result.begin(); dst != dstend; dst++) {
    *dst = scalar / *dst;
  }
  return result;
}

//----------------------jadp start of mods  13/7/2001

// gRowVector<bool, lhstortype> && gRowVector<bool, rhstortype>
template<typename lhstortype, typename rhstortype>
gRowVector<bool, DENSE>
operator&&(const gRowVector<bool, lhstortype>& vecl,
           const gRowVector<bool, rhstortype>& vecr) {

#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator&&(const gRowVector<bool, lhstortype>&, const gRowVector<bool, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  } // if
#endif // CHK_IDX

  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = ( vecl[indx] && vecr[indx] );
  } // while
  return vecret;
} // operator&&

// gRowVector<bool, lhstortype> || gRowVector<bool, rhstortype>
template<typename lhstortype, typename rhstortype>
gRowVector<bool, DENSE>
operator||(const gRowVector<bool, lhstortype>& vecl,
           const gRowVector<bool, rhstortype>& vecr) {

#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator||(const gRowVector<bool, lhstortype>&, const gRowVector<bool, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  } // if
#endif // CHK_IDX

  gRowVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] || vecr[indx]);
  } // while
  return vecret;
} // operator||

// gColVector<bool, lhstortype> && gColVector<bool, rhstortype>
template<typename lhstortype, typename rhstortype>
gColVector<bool, DENSE>
operator&&(const gColVector<bool, lhstortype>& vecl,
           const gColVector<bool, rhstortype>& vecr) {

#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator&&(const gColVector<bool, lhstortype>&, const gColVector<bool, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  } // if
#endif // CHK_IDX

  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] && vecr[indx]);
  } // while
  return vecret;
} // operator&&

// gColVector<bool, lhstortype> || gColVector<bool, rhstortype>
template<typename lhstortype, typename rhstortype>
gColVector<bool, DENSE>
operator||(const gColVector<bool, lhstortype>& vecl,
           const gColVector<bool, rhstortype>& vecr) {

#ifdef CHK_IDX
  if ( vecl.size() != vecr.size() ) {
    gErrorSize( "",
    "operator||(const gColVector<bool, lhstortype>&, const gColVector<bool, rhstortype>&)",
    "",
    vecl.size(), vecr.size() );
  } // if
#endif // CHK_IDX

  gColVector<bool, DENSE> vecret( vecl.size() );
  gLinearTypes::T_Index indx = vecl.size();
  while (indx--) {
    vecret[indx] = (vecl[indx] || vecr[indx]);
  } // while
  return vecret;
} // operator||

//----------------------jadp end of mods  13/7/2001

//----------------------lude start of mods  22/11/2001

// basic vector unary function
template<typename T, typename T_Stor, typename T_Functor>
typename gRowVector<T, T_Stor>::T_Dense
basicFunc(const gRowVector<T, T_Stor>& arg) {
  typedef typename gRowVector<T, T_Stor>::T_Dense T_Result;
  T_Result result( arg.size() );
  T_Functor functor;
  std::transform(arg.begin(), arg.end(), result.begin(), functor);
  return result;
}

// basic vector unary complex->real function
template<typename T, typename T_Stor, typename T_Functor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
complexToRealFunc(const gRowVector<T, T_Stor>& arg) {
  typedef typename gCppTraits<T>::T_Real T_Real;
  typedef gRowVector<T_Real, DENSE> T_Result;
  T_Result result( arg.size() );
  T_Functor functor;
  std::transform(arg.begin(), arg.end(), result.begin(), functor);
  return result;
}

// complex conjugate
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
conj(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcConj<T> >(arg);
}

// real component
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
real(const gRowVector<T, T_Stor>& arg) {
  return complexToRealFunc< T, T_Stor, funcReal<T> >(arg);
}

// imaginary component
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
imag(const gRowVector<T, T_Stor>& arg) {
  return complexToRealFunc< T, T_Stor, funcImag<T> >(arg);
}

// absolute value
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
abs(const gRowVector<T, T_Stor>& arg) {
  return complexToRealFunc< T, T_Stor, funcAbs<T> >(arg);
}

// argument (angle)
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
arg(const gRowVector<T, T_Stor>& arg) {
  return complexToRealFunc< T, T_Stor, funcArg<T> >(arg);
}

// norm (abs ^ 2)
template<typename T, typename T_Stor>
gRowVector<typename gCppTraits<T>::T_Real, DENSE>
norm(const gRowVector<T, T_Stor>& arg) {
  return complexToRealFunc< T, T_Stor, funcNorm<T> >(arg);
}

// cosine
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
cos(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcCos<T> >(arg);
}

// hyperbolic cosine
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
cosh(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcCosh<T> >(arg);
}

// sine
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sin(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcSin<T> >(arg);
}

// sinh
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sinh(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcSinh<T> >(arg);
}

// exp
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
exp(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcExp<T> >(arg);
}

// log
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
log(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcLog<T> >(arg);
}

// sqrt
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
sqrt(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcSqrt<T> >(arg);
}

// step function
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
step(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcStep<T> >(arg);
}

// impulse function
template<typename T, typename T_Stor>
typename gRowVector<T, T_Stor>::T_Dense
impulse(const gRowVector<T, T_Stor>& arg) {
  return basicFunc< T, T_Stor, funcImpulse<T> >(arg);
}

//----------------------lude end of mods  22/11/2001

// write gRowVector length, followed by contents to ostream
template<typename T, typename vecstortype>
std::ostream& operator<<(std::ostream& os,
       const gRowVector<T, vecstortype>& vec) {
  typedef typename gRowVector<T, vecstortype>::T_Index T_Index;
  // long oldFlags; JADP 1May2003
  std::ios::fmtflags oldFlags;
  T_Index N = vec.size();
  os << N << std::endl;
  oldFlags = os.flags(std::ios::scientific);
  for (T_Index indx = 0; indx < N; indx++) {
    os << vec[indx] << ' ';
  }
  os.flags(oldFlags);
  return os << std::endl;
}

// write gColVector length, followed by contents to ostream
template<typename T, typename vecstortype>
std::ostream& operator<<(std::ostream& os,
       const gColVector<T, vecstortype>& vec) {
  typedef typename gColVector<T, vecstortype>::T_Index T_Index;
  T_Index N = vec.size();
  os << N << std::endl;
  std::ios::fmtflags oldFlags = os.flags(std::ios::scientific);
  for (T_Index indx = 0; indx < N; indx++) {
    os << vec[indx] << ' ';
  }
  os.flags(oldFlags);
  return os << std::endl;
}

// read gRowVector length, followed by contents from istream
template<class T, typename vecstortype>
std::istream& operator>>(std::istream& is,
       gRowVector<T, vecstortype>& vec) {
  typedef typename gRowVector<T, vecstortype>::T_Index T_Index;
  T_Index N;
  is >> N;
  vec.resize(N);
  for (T_Index indx = 0; indx < N; indx++) {
    is >> vec[indx];
  }
  return is;
}

// read gColVector length, followed by contents from istream
template<typename T, typename vecstortype>
std::istream& operator>>(std::istream& is,
       gColVector<T, vecstortype>& vec) {
  typedef typename gColVector<T, vecstortype>::T_Index T_Index;
  T_Index N;
  is >> N;
  vec.resize(N);
  for (T_Index indx = 0; indx < N; indx++) {
    is >> vec[indx];
  }
  return is;
}

// write gRowVector length, followed by contents to ostream - binary
template<typename T, typename vecstortype>
std::ostream& global_write(std::ostream& os,
               const gRowVector<T, vecstortype>& vec) {
  typedef typename gRowVector<T, vecstortype>::T_Index T_Index;
  T_Index N = vec.size();
  os.write( reinterpret_cast<const char*>(&N), sizeof(N) );
  for (T_Index indx = 0; indx < N; indx++) {
    os.write( reinterpret_cast<const char*>(&vec[indx]), sizeof(T) );
  }
  return os;
}

// write gColVector length, followed by contents to ostream - binary
template<typename T, typename vecstortype>
std::ostream& global_write(std::ostream& os,
               const gColVector<T, vecstortype>& vec) {
  typedef typename gColVector<T, vecstortype>::T_Index T_Index;
  T_Index N = vec.size();
  os.write( reinterpret_cast<const char*>(&N), sizeof(N) );
  for (T_Index indx = 0; indx < N; indx++) {
    os.write( reinterpret_cast<const char*>(&vec[indx]), sizeof(T) );
  }
  return os;
}

// read gRowVector length, followed by contents from istream - binary
template<class T, typename vecstortype>
std::istream& global_read(std::istream& is,
              gRowVector<T, vecstortype>& vec) {
  typedef typename gRowVector<T, vecstortype>::T_Index T_Index;
  T_Index N;
  is.read( reinterpret_cast<char*>(&N), sizeof(N) );
  vec.resize(N);
  is.read( reinterpret_cast<char*>(vec.getStartAddr()), sizeof(T) * N );
  return is;
}

// read gRowVector length, followed by contents from istream - binary
template<class T, typename vecstortype>
std::istream& global_read(std::istream& is,
              gColVector<T, vecstortype>& vec) {
  typedef typename gColVector<T, vecstortype>::T_Index T_Index;
  T_Index N;
  is.read( reinterpret_cast<char*>(&N), sizeof(N) );
  vec.resize(N);
  is.read( reinterpret_cast<char*>(vec.getStartAddr()), sizeof(T) * N );
  return is;
}

//----------------------jadp start of mods  13/7/2001

// Forms a (possibly shorter) vector of the indices where gRowVector vec is true
template<typename T_Stor>
gRowVector<gLinearTypes::T_Index, DENSE> find(const gRowVector<bool, T_Stor>& vec) {

  // first count them
  gLinearTypes::T_Index cnt = 0;
  gLinearTypes::T_Index indx = vec.size();
  while (indx--) {
    if ( vec[indx] ) {
      cnt++;
    } // if
  } // while

  // form idx vector
  gRowVector<gLinearTypes::T_Index, DENSE> ret(cnt);
  indx = vec.size();
  cnt = ret.size();
  while (indx--) {
    if ( vec[indx] ) {
      ret[--cnt] = indx;
    } // if
  } // while

  return ret;
} // find

// Forms a (possibly shorter) vector of the indices where gColVector vec is true
template<typename T_Stor>
gColVector<gLinearTypes::T_Index, DENSE> find(const gColVector<bool, T_Stor>& vec) {

  // first count them
  gLinearTypes::T_Index cnt = 0;
  gLinearTypes::T_Index indx = vec.size();
  while (indx--) {
    if ( vec[indx] ) {
      cnt++;
    } // if
  } // while

  // form idx vector
  gColVector<gLinearTypes::T_Index, DENSE> ret(cnt);
  indx = vec.size();
  cnt = ret.size();
  while (indx--) {
    if ( vec[indx] ) {
      ret[--cnt] = indx;
    } // if
  } // while

  return ret;
} // find

//----------------------jadp end of mods  13/7/2001

GLINEAR_NAMESPACE_END
