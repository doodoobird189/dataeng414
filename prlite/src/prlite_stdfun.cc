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

 /*****************************************************************************
 *               |                                                            *
 * Author        | Dr J.A. du Preez, et. al.                                  *
 * Copyright     | University of Stellenbosch, all rights reserved.           *
 * Created       | 1995/09/14                                                 *
 * Description   | Various useful routines                                    *
 *               |                                                            *
 * Contact Info  |                                                            *
 *        [jadp] | Du Preez, Johan Dr                 <dupreez@dsp.sun.ac.za> *
 *        [goof] | Esterhuizen, Gerhard                  <goof@dsp.sun.ac.za> *
 *        [lude] | Schwardt, Ludwig                  <schwardt@ing.sun.ac.za> *
 *               |                                                            *
 *****************************************************************************/

// patrec headers
#include "prlite_stdfun.hpp"
#include "prlite_matrix.hpp"

// standard headers
#include <limits>  // numeric_limits
#include <cmath>  // ceil, frexp, ldexp, floor, log, sqrt

using namespace std;

namespace prlite{

/**************************** Constants & Defines      ***********************/

const float INFTY = numeric_limits<float>::infinity();

/************************** Type conversion (rounding) ***********************/

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
float fround(float x,
             int precision) {
  if ( x < -numeric_limits<int>::max() || x > numeric_limits<int>::max() ) {
    PRLITE_FAIL("round: " << x << " out of range");
  } // if
  float scale = 1;
  if (precision >= 0) {
    for (int i = 0; i < precision; i++) {
      scale *= 10.0f;
    } // for i
  } // if
  else {
    for (int i = 0; i > precision; i++) {
      scale /= 10.0f;
    } // for i
  } // else
  return static_cast<float>(static_cast<float>(round2long(x * scale)) / scale);
} // fround

// ----------------------------------------------------------------------------
int truncFloatToInt(float x) {
  if ( x < -numeric_limits<int>::max() || x > numeric_limits<int>::max() ) {
    PRLITE_FAIL("trunc: " << x << " out of range");
  } // if
  return static_cast<int>(x);
} // truncToFloat

/************************** Math functions             ***********************/

float cln(float x) {
  return x > 1E-30f ? log(x) : -999.9f;
} // cln

/*
Vector<float> log(const Vector<float>& arr) {
  Vector<float> result( arr.size() );
  for (int i = 0; i < arr.size(); i++) {
    result[i] = log( arr[i] );
  } // for i
  return result;
} // log
*/

// ----------------------------------------------------------------------------
int sign(float x) {
  if (x > 0) {
    return 1;
  } // if
  else if (x < 0) {
    return -1;
  } // else if
  return 0;
} // sign

// ----------------------------------------------------------------------------
float absfloat(float x) {
  return (x >= 0) ? x : -x;
} // absfloat

// ----------------------------------------------------------------------------
double absdouble(double x) {
  return (x >= 0) ? x : -x;
} // absdouble

// ----------------------------------------------------------------------------
int even(long val) {
  return ( (val / 2L) * 2L == val ) ? 1 : 0;
} // even


// ----------------------------------------------------------------------------
// projects x on the line connecting (x1, y1) and (x2, y2)
float linInterpolate(float x1,
                     float y1,
                     float x2,
                     float y2,
                     float x) {
  return y1 + (x - x1) / (x2 - x1) * (y2 - y1);
} // linInterpolate

// ----------------------------------------------------------------------------
// Least-squares polynomial curve fitting   (RvdM 18/09/97)
// x = independant variable sequence
// y = dependant variable sequence
// m = order of polynomial
// coef = polynomial coefficients [a_0, a_1, ... , a_m] (return value)
// y2 = evaluation of polynomial at x values
// returns the polynomial regression error
float polyFit(const gLinear::gRowVector<float>& x,
        const gLinear::gRowVector<float>& y,
        int m,
        gLinear::gRowVector<float>& coef) {

  int n = x.size();

  if (n <= m) {
    PRLITE_FAIL("(polyFit) Error ! Data matrix underdetermined. Number of data points" << endl
         << "                  less than or equal to order of polynomial.");
  } // if

  gLinear::gRowVector<float> b(m + 1);
  gLinear::gRowMatrix<float> V(m + 1, m + 1);
  gLinear::gRowVector<float> t(2 * m + 1);
  gLinear::gRowVector<float> xc(n);
  float dummy, xd;
  float error = 0.0;

  xc.assignToAll(1.0);
  t.assignToAll(0.0);
  b.assignToAll(0.0);
  t[0] = static_cast<float>(n);

  // setup V and b matrices
  for (int i = 1; i < 2 * m + 1; i++) {
    for (int j = 0; j < n; j++) {
      if ( (i - 1) < m + 1 ) {
        b[i - 1] += y[j] * xc[j];
      } // if
      xc[j] *= x[j];
      t[i] += xc[j];
    } // for
  } // for
  for (int i = 0; i < m + 1; i++) {
    for (int j = 0; j < m + 1; j++) {
      V[i][j] = t[i + j];
    } // for
  } // for

  int fail;
  gLinear::gRowMatrix<float> iV = inv(V, fail);
  PRLITE_ASSERT(!fail, "Failed to invert matrix.");

  coef = b * iV;
  gLinear::gRowVector<float> y2(n);

  // evaluate polynomial
  for (int i = 0; i < n; i++) {
    y2[i] = coef[0];
    xd = 1.0;
    for (int j = 1; j < m + 1; j++) {
      xd = xd * x[i];
      y2[i] = y2[i] + coef[j] * xd;
    } // for
    dummy = y[i] - y2[i];
    error = error + dummy * dummy;
  } // for

  error = sqrt( error / static_cast<float>( n - (m + 1) ) );

  return error;

} // polyFit

// ----------------------------------------------------------------------------
gLinear::gRowVector<float> polyEval(const gLinear::gRowVector<float>& x,
                       const gLinear::gRowVector<float>& coef) {

  int n = x.size();
  int m = coef.size() - 1;
  gLinear::gRowVector<float> y(n);
  float xd;

  // evaluate polynomial
  for (int i = 0; i < n; i++) {
    y[i] = coef[0];
    xd = 1.0;
    for (int j = 1; j < m + 1; j++) {
      xd = xd * x[i];
      y[i] = y[i] + coef[j] * xd;
    } // for
    if ( !(y[i] <= INFTY) ) {
      cout << "PROBLEM TIME AGAIN !!!! \n\n";
      cout << "y = " << y[i];
      cout << "x = " << x[i] << endl;
      cout << "coef = " << coef << endl;
    } // if
  } // for
  return y;
} // polyEval

// ----------------------------------------------------------------------------

float sortedLookup( float x,
                    const gLinear::gRowVector<float>& sortedX,
                    const gLinear::gRowVector<float>& outVals,
                    unsigned interp) {

  unsigned i = 0;
  unsigned j = sortedX.size() -1;
  unsigned t = (i + j)/2;

  // we intentionally do no extrapolation to work properly with
  // probabilistic lookups like PAV
  if (x <= sortedX[0]) { return outVals[0]; } // if
  if (x >= sortedX[j]) { return outVals[j]; } // if

  do { // do a binary search to find the position of x
    if ( x < sortedX[t] ){ // x is lower down
      j = t;
      t = (i + j)/2;
      if (i == t) break;
    } // if
    else if (x > sortedX[t]) { // x is further up
      i = t;
      t = (i + j + 1)/2;
      if (j == t) break;
    } // else if
    else { // precise match
      return outVals[t];
    } // else
  } while (1); // while

  if (interp) { // somewhere between sortedX[i] and sortedX[j]
    return linInterpolate(sortedX[i], outVals[i], sortedX[j], outVals[j], x);
  } // if
  else {
    return outVals[i];
  } // else
} // sortedLookup


// ----------------------------------------------------------------------------
int nextPowerOf2(unsigned from) {
  int exp2;
  double frac2;
  frac2 = frexp(static_cast<double>(from), &exp2);
  if (frac2 == 0.5) {
    exp2--;
  } // if
  return round2int( ldexp(1.0, exp2) );
} // nextPowerOf2

// ----------------------------------------------------------------------------
int intPow(const int base,
           const unsigned int exponent) {
  int ans(1);
  for (unsigned int i(1); i <= exponent; i++) {
    ans *= base;
  } // for i
  return ans;
} // intPow

/************************** Frame indexing             ***********************/

/*!!! Fout

  t_mid(0) is die middelpunt van raam 0. tS is die tydspasiering tussen rame.
  t_mid(i) = t_mid(0) + tS * i

  ------------------------ t_old -----------------
  Neem i_old(t) as die huidige def vir die raam indeks by tyd t. Verder is vas
  gekies t_mid(0) = tS:

  i_old(t) = floor( ( t - t_mid(0) ) / tS )
           = floor( (t - tS) / tS )

  dan volg die huidige def vir die tyd van raam i:

  t_old(i) = t_mid(0) + tS * i met vaste aanname dat
           = tS * (i + 1)

  Vir die grens tusen twee rame sien dat
  i_old( t_mid(i) ) = i_old(t_mid(0) + tS * i)
                         = floor( (t_mid(0) + tS * i - tS) / tS )
                         = floor(tS * i / tS)   as t_mid(0) == tS
                         = floor(i)
                         = i, terwyl
  i_mid(t_mid(i)-) = floor(i-)
                   = i - 1.
  Met ander woorde die raam domineer die area (let op die oop en geslote grense)
  [ t_mid(i) .. t_mid(i) + tS ), oftewel die tweede helfte van die werklike raam

  ------------------------ t_mid -----------------
  Neem i_mid(t) as die raam waarvan die middelpunt die naaste aan die tyd t is.
  As dit presies in die middel van twee middelpunte is, vat ons die tweede een:

  i_mid(t) = round( ( t - t_mid(0) ) / tS )

  dan is die tyd in die middel van raam i:

  t_mid(i) = t_mid(0) + tS * i

  Vir die grens tusen twee rame sien dat
  i_mid(t_mid(i) + tS/2) = i_mid(t_mid(0) + tS * i + tS/2)
                       = round( ( t_mid(0) + tS * i + tS/2 - t_mid(0) ) / tS )
                       = round(tS * (i + 0.5) / tS)
                       = round(i + 0.5)
                       = i + 1, terwyl
  i_mid(t_mid(i) + tS/2-) = round(i + 0.5-)
                         = i.
  Met ander woorde die raam domineer die area (let op die oop en geslote grense)
  [ t_mid(i) - tS/2 .. t_mid(i) + tS/2 ).

  ------------------------ t_beg -----------------

  As ons in terme van die begintyd van 'n raam dink, dan is die punt halfpad
  tussen twee rame se middelpunte die begintyd van die tweede een van die twee.
  Voor hierdie punt is die eerste raam dominant. Met hierdie oogpunt dan is die
  raam op tyd t:

  i_beg(t) = floor( (t - t_mid(0) + tS/2) / tS )

  en die begintyd van hierdie raam:

  t_beg(i) = t_mid(i) - tS/2 = t_mid(0) + tS(i - 0.5)

  Vir die grens tusen twee rame i.t.v hierdie tyd, sien dat
  i_beg( t_beg(i) ) = i_beg(t_mid(i) - tS/2)
                    = i_beg(t_mid(0) + tS * i - tS/2)
                    = floor( (t_mid(0) + tS * i - tS/2 - t_mid(0) + tS/2) / tS )
                    = floor(tS * i / tS)
                    = floor(i)
                    = i, terwyl
  i_beg(t_beg(i)-) = floor(i+)
                    = i - 1  ???

  Met ander woorde die raam domineer die area van
  [ t_mid(i) - tS/2 .. t_mid(i) + tS/2).
  Lyk dus selfde as die t_mid definisie, maar raps goedkoper

  ------------------------ t_end -----------------

  As ons in terme van die eindtyd van 'n raam dink, dan is die punt halfpad tussen
  twee opeenvolgende rame se middelpunte die eindtyd vir die eerste van die twee.
  Verby hierdie punt is die tweede raam dominant. Met hierdie oogpunt dan is die
  raam op tyd t:

  i_end(t) = ceil( (t - t_mid(0) - tS/2) / tS )

  en die eindtyd van hierdie raam:

  t_end(i) = t_mid(i) + tS/2 = t_mid(0) + tS(i + 0.5)

  Vir die grens tusen twee rame i.t.v hierdie tyd, sien dat
  i_end( t_end(i) ) = i_end(t_mid(i) + tS/2)
                    = i_end(t_mid(0) + tS * i + tS/2)
                    = ceil( (t_mid(0) + tS * i + tS/2 - t_mid(0) - tS/2) / tS )
                    = ceil(tS * i / tS)
                    = ceil(i)
                    = i, terwyl
  i_end(t_end(i)+) = ceil(i+)
                    = i + 1
  Met ander woorde die raam domineer die area van
  (t_mid(i) - tS/2 .. t_mid(i) + tS/2].
  Dus, vir eindtyd gebaseerde transkripsie lyk t_end(i) na 'n meer
  korrekte definisie.

  ------------------------------------------------

  Let op dat
  i_beg(t_mid(i) + tS/2) = i + 1
  terwyl
  i_end(t_mid(i) + tS/2) = i
  wat reg lyk. Die raam wat op tyd t_mid(i) + tS/2 eindig is die huidige raam i
  terwyl die een wat daar begin is die volgende raam i + 1.

  Onderstaande formules is net
  reg met t_mid(0) == frameShift - normaalweg die geval.

*/

/*
  i_mid(t) = round( ( t - t_mid(0) ) / tS )
  t_mid(i) = t_mid(0) + tS * i
*/

int i_mid(double t, double tS, double tm0) {
  return round2int( (t - tm0) / tS + 1e-8 );
} // i_mid

double t_mid(int i, double tS, double tm0) {
  return tm0 + tS * i;
} // t_mid

void printmid(double t, double tS, double tm0) {
  cout << "\tm: "
       << "\t " << i_mid(t, tS, tm0)
       << "\t " << t_mid(i_mid(t, tS, tm0), tS, tm0)
       << "\t " << i_mid(t_mid(i_mid(t, tS, tm0), tS, tm0), tS, tm0);
} // printmid

/*
  i_beg(t) = floor( (t - t_mid(0) + tS/2) / tS )
  t_beg(i) = t_mid(i) - tS/2 = t_mid(0) + tS(i - 0.5)
*/

int i_beg(double t, double tS, double tm0) {
  return int( floor( ( t - tm0 + tS/2) / tS + 1e-8 ) );
} // i_beg

double t_beg(int i, double tS, double tm0) {
  return tm0 + tS * (i - 0.5);
} // t_beg

void printbeg(double t, double tS, double tm0) {
  cout << "\tb: "
       << "\t " << i_beg(t, tS, tm0)
       << "\t " << t_beg(i_mid(t, tS, tm0), tS, tm0)
       << "\t " << i_beg(t_mid(i_mid(t, tS, tm0), tS, tm0), tS, tm0);
} // printbeg

/*
  i_end(t) = ceil( (t - t_mid(0) - tS/2) / tS )
  t_end(i) = t_mid(i) + tS/2 = t_mid(0) + tS(i + 0.5)
*/

int i_end(double t, double tS, double tm0) {
  return int( ceil( ( t - tm0 - tS/2) / tS + 1e-8 ) );
} // i_end

double t_end(int i, double tS, double tm0) {
  return tm0 + tS * (i + 0.5);
} // t_end

void printend(double t, double tS, double tm0) {
  cout << "\te: "
       << "\t " << i_end(t, tS, tm0)
       << "\t " << t_end(i_mid(t, tS, tm0), tS, tm0)
       << "\t " << i_end(t_mid(i_mid(t, tS, tm0), tS, tm0), tS, tm0);
} // printend

/*
  // toets vir t_beg, t_mid en t_end defs
    double tm0 = 10.0;
    double tS = 10.0;
    double eps = 0.1;
    for (double t = 0; t < 100; t+= tS/2.0) {
      cout << t-eps << ":"; printbeg(t-eps, tS, tm0); printmid(t-eps, tS, tm0); printend(t-eps, tS, tm0); cout << endl;
      cout << t+0.0 << ":"; printbeg(t+0.0, tS, tm0); printmid(t+0.0, tS, tm0); printend(t+0.0, tS, tm0); cout << endl;
      cout << t+eps << ":"; printbeg(t+eps, tS, tm0); printmid(t+eps, tS, tm0); printend(t+eps, tS, tm0); cout << endl;
    } // for t
 */

//  i_end(t) = ceil( (t - t_mid(0) - tS/2) / tS )
int frameIdx(double t,
             double tS) {
  double t_mid0 = tS; // for the moment assume t_mid(0) == tS
  return i_mid(t, tS, t_mid0);
  // maybe a bit pedantic, but within accuracy the right hand boundary is included.
  //return static_cast<int>( ceil( (t - t_mid0 - tS/2) / tS -1E-8 ) ); // !!!printend: die -1E-8 onstabiel op linkergrens
} // frameIdx

// t_end(i) = t_mid(i) + tS/2 = t_mid(0) + tS(i + 0.5)
double frameTime(int i,
                 double tS) {
  double t_mid0 = tS; // for the moment assume t_mid(0) == tS
  return t_mid(i, tS, t_mid0);
  //return t_mid0 + tS * (0.5 + i);
} // frameTime

/* Hierdie is die ou weergawe - sien i_old beskrywing hierbo
 * Aanvaar eerste raam se middelpunt is by frameShift en opeenvolgende rame is
 * met intervalle frameShift gespasieer.
 * frameIdx: verskaf die naaste raam met middelpunt links-van/ op  time, dis meer sinvol
 *           vir eindtye aangesien daar maksimaal 'n halwe raam se data van anderkant
 *           die grens kan kom. Onthou, dis eindtye, etime = 0 sal raamno -1 of selfs
 *           vroeer lewer.
 * frameTime: lewer die tyd by die middel van die raam
 */
// Changed floats in frameIdx to doubles due to precision problems. (RvdM - 29/05/97)
int frameIdxOld(float etime,
                float frameShift) {
  return static_cast<int>( floor( ( static_cast<double>(etime) - static_cast<double>(frameShift) ) / static_cast<double>(frameShift) ) );
} // frameIdxOld

// ----------------------------------------------------------------------------
float frameTimeOld(int idx,
                   float frameShift) { // tyd in middel van raam
  return static_cast<float>(frameShift * static_cast<float>(idx) + frameShift);
} // frameTimeOld

// ----------------------------------------------------------------------------
unsigned noOfVectorsInFrames(unsigned noOfFrames,
                             unsigned frmLen,
                             unsigned frmShift) {
  if ( noOfFrames == numeric_limits<unsigned int>::max() ) {
    return numeric_limits<unsigned int>::max();
  } // if
  if (noOfFrames == 0) {
    return 0;
  } // if
  return frmLen + (noOfFrames - 1) * frmShift;
} // noOfVectorsInFrames

} // prlite

/*
// naaste raam met middelpunt links van time; beter vir eindtye
// onthou, dis eindtye, etime = 0 sal raamno -1 of selfs vroeer lewer.
int frameIdx(float etime,
             float frameDuration,
             float frameShift) {
  // return round( (time - frameDuration / 2.0) / frameShift );
  return floor( (etime - frameDuration / 2.0) / frameShift );
} // frameIdx

float frameTime(int idx,
                float frameDuration,
                float frameShift) {
  return(frameShift * idx + frameDuration / 2.0);
} // frameTime
*/
