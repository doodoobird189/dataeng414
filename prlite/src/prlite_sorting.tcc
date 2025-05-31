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
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : pre-2000
 * Copyright  : University of Stellenbosch, all rights retained
 */

// standard headers
#include <cmath>  // ceil, floor
#include <algorithm>  // max, min
#include <iostream>
#include <random>

namespace prlite{

// Quicksort: Uncomment for median pivot, otherwise a random pivot is used
// #define QUICKSORT_MEDIAN_PIVOT

template<typename VEC>
void quickSortDown(const VEC& vec,
                   std::vector<int>& index,
                   int from,
                   int upto) {
  if (upto < 0) {
    upto = vec.size() - 1;
  } // if
  int i(from);
  int j(upto);
  int swapi;

#ifndef QUICKSORT_MEDIAN_PIVOT
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rgenInt(from, upto);
  typename VEC::value_type pivot( vec[ index[ rgenInt(gen) ] ] );
  //typename VEC::value_type pivot( vec[ index[ rgenInt( emdw::randomEngine() ) ] ] );
#else
  typename VEC::value_type pivot( vec[ index[ (from + upto) / 2 ] ] );
#endif // QUICKSORT_MEDIAN_PIVOT

  do {
    while (vec[ index[i] ] > pivot) {
      i++;
    } // while
    while ( pivot > vec[ index[j] ] ) {
      j--;
    } // while
    if (i <= j) {
      swapi = index[i];
      index[i++] = index[j];
      index[j--] = swapi;
    } // if
  } while (i <= j);
  if (from < j) {
    quickSortDown(vec, index, from, j);
  } // if
  if (i < upto) {
    quickSortDown(vec, index, i, upto);
  } // if
} // quickSortDown


template<typename VEC>
std::vector<int> sortDown(const VEC& vec,
                     int from,
                     int upto) {
  typename VEC::size_type sz( vec.size() );
  std::vector<int> sorted(sz);
  if (upto < 0) {
    upto = vec.size() - 1;
  } // if
  for (int i = from; i <= upto; i++) {
    sorted[i] = i;
  } // for i

  quickSortDown(vec, sorted, from, upto);
  return sorted;
} // sortDown


template<typename VEC>
void quickSortUp(const VEC& vec,
                 std::vector<int>& index,
                 int from,
                 int upto) {
  if (upto < 0) {
    upto = vec.size() - 1;
  } // if
  int i(from);
  int j(upto);
  int swapi;

#ifndef QUICKSORT_MEDIAN_PIVOT
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> rgenInt(from, upto);
  typename VEC::value_type pivot( vec[ index[ rgenInt(gen) ] ] );
  //typename VEC::value_type pivot( vec[ index[ rgenInt( emdw::randomEngine() ) ] ] );
#else
  typename VEC::value_type pivot( vec[ index[ (from + upto) / 2 ] ] );
#endif // QUICKSORT_MEDIAN_PIVOT

  do {
    while (vec[ index[i] ] < pivot) {
      i++;
    } // while
    while ( pivot < vec[ index[j] ] ) {
      j--;
    } // while
    if (i <= j) {
      swapi = index[i];
      index[i++] = index[j];
      index[j--] = swapi;
    } // if
  } while (i <= j);
  if (from < j) {
    quickSortUp(vec, index, from, j);
  } // if
  if (i < upto) {
    quickSortUp(vec, index, i, upto);
  } // if
} // quickSortUp


template<typename VEC>
std::vector<int> sortUp(const VEC& vec,
                   int from,
                   int upto) {
  int sz( vec.size() );
  std::vector<int> sorted(sz);
  if (upto < 0) {
    upto = vec.size() - 1;
  } // if
  for (int i = from; i <= upto; i++) {
    sorted[i] = i;
  } // for i

  quickSortUp(vec, sorted, from, upto);
  return sorted;
} // sortUp


template<typename VEC>
void partialSort(const VEC& vec,
                    const std::vector<int>& toFind,
                    std::vector<int>& index,
                    int from,
                    int upto) {
  if (upto < 0) {
    upto = vec.size() - 1;
  } // if
  int i(from);
  int j(upto);

  typename VEC::value_type median( vec[ index[ (from + upto) / 2 ] ] );
  int swapi;

  do { // order elements around median
    while (vec[ index[i] ] < median) { // skip the leftside ones already smaller
      i++;
    } // while
    while ( median < vec[ index[j] ] ) { // skip the rightside ones already bigger
      j--;
    } // while
    if (i <= j) { // swop the two out of place values
      swapi = index[i];
      index[i++] = index[j];
      index[j--] = swapi;
    } // if
  } while (i <= j);

  // use a quicksort, but only on the appropriate segments
  for (unsigned int k = 0; k < toFind.size(); k++) {
    if (from < j && from <= toFind[k] && toFind[k] <= j) {
      partialSort(vec, toFind, index, from, j);
      break;
    } // if
  } // for k
  for (unsigned int k = 0; k < toFind.size(); k++) {
    if (i < upto && i <= toFind[k] && toFind[k] <= upto) {
      partialSort(vec, toFind, index, i, upto);
      break;
    } // if
  } // for k

} // partialSort


// the old one
template<typename T>
T percentile_old(const std::vector<T>& vec,
             float p,
             int from,
             int upto) { // p as a fraction

  // Kreyszig bl 1244
  if (upto <= 0) {
    upto = vec.size() - 1;
  } // if
  int len(upto - from + 1);

  std::vector<int> toFind(2);
  double actualIdx = from + p * static_cast<double>(len - 1);
  double lowerIdx = std::max(std::floor(actualIdx), 0.0);
  double upperIdx = std::min(std::ceil(actualIdx), static_cast<double>(upto) );
  //std::cout << "Lower " << lowerIdx << " Actual " << actualIdx << " Upper " << upperIdx << std::endl;
  toFind[0] = static_cast<int>(lowerIdx);
  toFind[1] = static_cast<int>(upperIdx);

  std::vector<int> idx( vec.size() );
  for (int i = from; i <= upto; i++) {
    idx[i] = i;
  } // for
  partialSort(vec, toFind, idx, from, upto);
  //std::cout << "lower val: " << vec[ idx[ toFind[0] ] ] << " and upper val: " << vec[ idx[ toFind[1] ] ] << std::endl;
  return static_cast<T>(
  toFind[0] == toFind[1] ?
    vec[ idx[ toFind[0] ] ] :
    ( vec[ idx[ toFind[0] ] ] *( upperIdx - actualIdx) +
           vec[ idx[ toFind[1] ] ] * (actualIdx - lowerIdx ) )
           );
} // percentile_old

// this one is the correct one for future use. main difference is to
// realise that with N points, the correct percentile value of the
// n'th point is (n-0.5)/N, i.e. half of that n'th point goes to the
// lower interval and half to the higher interval. Between two points
// we split the percentile value linearly, i.e. we will always have a
// unique answer for percentile values that fall between two actual
// data points. At the endpoints we limit the values to 0 and 1,
// i.e. we do not extrapolate beyond the endpoints.

template<typename VEC>
typename VEC::value_type percentile(const VEC& vec, double r, int
from, int upto) { // p as a fraction

  // Kreyszig bl 1244
  if (upto <= 0) {
    upto = vec.size() - 1;
  } // if
  int len(upto - from + 1);

  std::vector<int> toFind(2);
  double actualIdx = from + r * static_cast<double>(len)-0.5;
  double floorIdx = std::max(std::floor(actualIdx), 0.0);
  double ceilIdx = std::min(std::ceil(actualIdx), static_cast<double>(upto) );
  //std::cout << "Lower " << floorIdx << " Actual " << actualIdx << " Upper " << ceilIdx << std::endl;
  toFind[0] = static_cast<int>(floorIdx);
  toFind[1] = static_cast<int>(ceilIdx);

  std::vector<int> idx( vec.size() );
  for (int i = from; i <= upto; i++) {
    idx[i] = i;
  } // for
  partialSort(vec, toFind, idx, from, upto);
  //std::cout << "lower val: " << vec[ idx[ toFind[0] ] ] << " and upper val: " << vec[ idx[ toFind[1] ] ] << std::endl;

  // !!! c++11 dink die float moet VEC::value_type wees, maar compile nie
  return static_cast<float>(toFind[0] == toFind[1] ?
    vec[ idx[ toFind[0] ] ] :
    ( vec[ idx[ toFind[0] ] ] *( ceilIdx - actualIdx) +
      vec[ idx[ toFind[1] ] ] * (actualIdx - floorIdx ) ));
} // percentile_new


template<typename VEC>
void  percentileRange(const VEC& vec,
                      double rLower,
                      double rUpper,
                      typename VEC::value_type& percLower,
                      typename VEC::value_type& percUpper,
                      int from,
                      int upto) {
  if (upto <= 0) {
    upto = vec.size() - 1;
  } // if
  int len(upto - from + 1);

  std::vector<int> toFind(4);

  // first the lower
  double lowerActualIdx = from + rLower * static_cast<double>(len)-0.5;
  double lowerFloorIdx = std::max(std::floor(lowerActualIdx), 0.0);
  double lowerCeilIdx = std::min(std::ceil(lowerActualIdx), static_cast<double>(upto) );
  // cout << "Lower " << floorIdx << " Actual " << actualIdx << " Upper " << ceilIdx << endl;
  toFind[0] = static_cast<int>(lowerFloorIdx);
  toFind[1] = static_cast<int>(lowerCeilIdx);

  // then the upper
  double upperActualIdx = from + rUpper * static_cast<double>(len)-0.5;
  double upperFloorIdx = std::max(std::floor(upperActualIdx), 0.0);
  double upperCeilIdx = std::min(std::ceil(upperActualIdx), static_cast<double>(upto) );
  // cout << "Lower " << floorIdx << " Actual " << actualIdx << " Upper " << ceilIdx << endl;
  toFind[2] = static_cast<int>(upperFloorIdx);
  toFind[3] = static_cast<int>(upperCeilIdx);

  // now do the searching
  std::vector<int> idx( vec.size() );
  for (int i = from; i <= upto; i++) {
    idx[i] = i;
  } // for
  partialSort(vec, toFind, idx, from, upto);

  percLower = toFind[0] == toFind[1] ?
    vec[ idx[ toFind[0] ] ] :
    ( vec[ idx[ toFind[0] ] ] *( lowerCeilIdx - lowerActualIdx) +
      vec[ idx[ toFind[1] ] ] * (lowerActualIdx - lowerFloorIdx ) );

  percUpper = toFind[2] == toFind[3] ?
    vec[ idx[ toFind[2] ] ] :
    ( vec[ idx[ toFind[2] ] ] *( upperCeilIdx - upperActualIdx) +
      vec[ idx[ toFind[3] ] ] * (upperActualIdx - upperFloorIdx ) );

} // percentileRangeNew

template<typename VEC>
typename VEC::value_type  median(const VEC& vec,
                                 int from,
                                 int upto) {
  return percentile(vec, 0.5, from, upto);
} // median

template<typename VEC>
VEC medianFilter(const VEC& vec,
                 int width) {
  VEC output( vec.size() );
  int last(vec.size() - 1);
  int firstFull(width / 2);
  int lastFull(last - firstFull);

  // first do the bits at the ends which must use lower order filters.
  output[0] = vec[0];
  output[last] = vec[last];
  for (int i = 1; i < firstFull; i++) {
    output[i] = median(vec, 0, 2 * i);
    output[last - i] = median(vec, last - 2 * i, last);
  } // for i

  // now for the ones that use the full filter
  for (int i = firstFull; i <= lastFull; i++) {
    output[i] = median(vec, i - firstFull, i + firstFull);
  } // for i

  return output;
} // medianFilter

} // namespace prlite
