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
  DEPRECATED! rather use the new C++17 functions
*/

/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */


#ifndef PRLITE_SORTING_HPP
#define PRLITE_SORTING_HPP

// patrec headers
#include <vector>

namespace prlite{

/**
 * Comparison function for use with STL sort and stable_sort, and also
 * STL containers.

How to use:
<PRE>
  const vector<string>& stuffToSort = ...;
  ...

  vector<int> indices( stuffToSort.size() );
  for (unsigned ii=0; ii<stuffToSort.size(); ii++)
      indices[ii] = ii;

  std::sort(stuffToSort.begin(), stuffToSort.end(),
            IndexCompare<vector<string> >(stuffToSort));
</PRE>
  This also works with gLinear Vector's iterators.  stuffToSort can be any container that supports []
 */
template<typename T>
struct IndexCompare {
  const T& storage;
  IndexCompare(const T& strg) : storage(strg) { }
  bool operator()(typename T::size_type i1, typename T::size_type i2) {
    return storage[i1] < storage[i2];
  }
};


/** quickSortDown:
  vec        : the vector to sort
  index      : Input: same length as vec and contains indices into vec. Initially
               index[from]..index[upto] must contain values from, from + 1, ... upto.
               Return: index[from]..index[upto] contains indices for a sorted vec.
               Outside this range it is undefined.
  from, upto : range in sorted over which sorting will be done. upto < 0 means
               upto = vec.size() - 1
*/

template<typename VEC>
void quickSortDown(const VEC& vec,
                   std::vector<int>& index,
                   int from = 0,
                   int upto = -1);

/** sortDown:
  vec        : the vector to sort
  from, upto : range in sorted over which sorting will be done. upto < 0 means
               upto = vec.size() - 1
*/

template<typename VEC>
std::vector<int> sortDown(const VEC& vec,
                     int from = 0,
                     int upto = -1);

/** quickSortUp:
  vec        : the vector to sort
  index      : Input: same length as vec and contains indices into vec. Initially
               index[from]..index[upto] must contain values from, from + 1, ... upto.
               Return: index[from]..index[upto] contains indices for a sorted vec
  from, upto : range in sorted over which sorting will be done. upto < 0 means
               upto = vec.size() - 1
*/

template<typename VEC>
void quickSortUp(const VEC& vec,
                 std::vector<int>& index,
                 int from = 0,
                 int upto = -1);

/** sortUp:
  vec        : the vector to sort
  from, upto : range in sorted over which sorting will be done. upto < 0 means
               upto = vec.size() - 1
*/

template<typename VEC>
std::vector<int> sortUp(const VEC& vec,
                   int from = 0,
                   int upto = -1);

/**
 * @brief performs a partial sort over vec[from:upto].
 *
 * @param vec the container under investigation. must support [i] and
 *   size() and the value_type trait
 *
 * @param toFind this routine can determine a number of N'th smallest
 *   values simultaneously. toFind lists them. NOTE: N = from
 *   specifies the smallest while N = upto specifies the biggest.
 * Example: if you want to search vec[20:40] for 2nd smallest and 3rd
 *   largest values, you will set index = [21,38]
 *
 * @param index On entry: same length as vec and contains indices into
 *   vec. Initially index[from:upto]  = [from:1:upto]
 *
 * @param from:upto Range of indices in vec over which partial sorting
 *   will be done. Sections of vec outside this range is left untouched.
 *   upto = -1 is taken to mean upto the last one.
 *
 * @return by reference an "index" vector of the same length as vec
 * such that:
 *
 * Unsorted ranges:
 *   index[0:from-1] remains untouched, i.e. contain indices 0:from-1.
 *   index[upto+1:last] remain untouched, i.e. contain indices upto+1:last
 * Partially sorted ranges:
 *   vec[from:index[toFind[j]-1]  <= (leq) vec[index[toFind[j]] and
 *   vec[index[toFind[j]+1]:upto] >= (geq) vec[index[toFind[j]].
 *
 */
template<typename VEC>
void partialSort(const VEC& vec,
                    const std::vector<int>& toFind,
                    std::vector<int>& index,
                    int from = 0,
                    int upto = -1);

/// the old one
template<typename T>
T percentile_old(const std::vector<T>& vec,
             float p,
             int from = 0,
             int upto = -1);

/**
 * @brief find the value of the requested percentile of a sequence of values
 *
 * @param vec the container under investigation. must support [i] and
 * size() and the value_type trait
 *
 * @param r The ratio to calculate the percentile at - expressed as a
 * FRACTION (i.e. not percentage)
 *
 * @param from First index in vec range over which search takes place
 * @param upto Last index in vec range over which search takes place
 */
template<typename VEC>
typename VEC::value_type  percentile(const VEC& vec,
                                     double r,
                                     int from = 0,
                                     int upto = -1);

/**
 * @brief find a lower and upper percentile point simultaneously. This
 * is slightly more efficient than calling percentile twice.
 *
 * @param vec the container under investigation. must support [i] and
 * size() and the value_type trait
 *
 * @param rLower The ratio to calculate the lower percentile at - expressed as a FRACTION
 * (i.e. not percentage)
 *
 * @param rUpper The ratio to calculate the upper percentile at - expressed as a FRACTION
 * (i.e. not percentage)
 *
 * @param percLower Returns the lower percentile
 *
 * @param percUpper Returns the upper percentile
 *
 * @param from First index in vec range over which search takes place
 * @param upto Last index in vec range over which search takes place
 */
template<typename VEC>
void  percentileRange(const VEC& vec,
                      double rLower,
                      double rUpper,
                      typename VEC::value_type& percLower,
                      typename VEC::value_type& percUpper,
                      int from = 0,
                      int upto = -1);

/**
 * @brief find the median of a sequence of values
 *
 * @param vec the container under investigation. must support [i] and
 * size() and the value_type trait
 *
 * @param from First index in vec range over which search takes place
 * @param upto Last index in vec range over which search takes place
 */
template<typename VEC>
typename VEC::value_type  median(const VEC& vec,
                                 int from = 0,
                                 int upto = -1);

/// a non-causal median filter - no lag; width must be odd
template<typename VEC>
VEC medianFilter(const VEC& vec,
                 int width);

} // namespace prlite

#include "prlite_sorting.tcc"

#endif // PRLITE_SORTING_HPP
