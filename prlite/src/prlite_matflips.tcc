/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2013
          PURPOSE: routines to flip matrices and vectors in various ways
          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

namespace prlite{

template <typename AssignType>
gLinear::gRowVector<AssignType>
flip(
  const gLinear::gRowVector<AssignType>& inVec){
  int len = inVec.size();
  gLinear::gRowVector<AssignType> outVec(len);
  for (unsigned i = 0; i < len; i++) {
    outVec[i] = inVec[len-1-i];
  } // for i
  return outVec;
} // flip

template <typename AssignType>
gLinear::gRowMatrix<AssignType>
rot90(
  const gLinear::gRowMatrix<AssignType>& inMat) {

  unsigned nRows = inMat.cols(); // new rows
  unsigned nCols = inMat.rows(); // new cols
  gLinear::gRowMatrix<AssignType> outMat(nRows, nCols);
  for (unsigned r = 0; r < nRows; r++) { // new indices
    for (unsigned c = 0; c < nCols; c++) {
      outMat(r,c) = inMat(nCols-1-c, r);
    } // for
  } // for
  return outMat;
} // rot90

template <typename AssignType>
gLinear::gRowMatrix<AssignType>
fliplr(
  const gLinear::gRowMatrix<AssignType>& inMat) {

  unsigned nRows = inMat.rows(); // new rows
  unsigned nCols = inMat.cols(); // new cols
  gLinear::gRowMatrix<AssignType> outMat(nRows, nCols);
  for (unsigned r = 0; r < nRows; r++) { // new indices
    for (unsigned c = 0; c < nCols; c++) {
      outMat(r,c) = inMat(r, nCols-1-c);
    } // for
  } // for
  return outMat;
} // fliplr

template <typename AssignType>
gLinear::gRowMatrix<AssignType>
flipud(
  const gLinear::gRowMatrix<AssignType>& inMat) {

  unsigned nRows = inMat.rows(); // new rows
  unsigned nCols = inMat.cols(); // new cols
  gLinear::gRowMatrix<AssignType> outMat(nRows, nCols);
  for (unsigned r = 0; r < nRows; r++) { // new indices
    for (unsigned c = 0; c < nCols; c++) {
      outMat(r,c) = inMat(nRows-1-r, c);
    } // for
  } // for
  return outMat;
} // flipud

template <typename AssignType>
gLinear::gRowMatrix<AssignType>
fliplrud(
  const gLinear::gRowMatrix<AssignType>& inMat) {

  unsigned nRows = inMat.rows(); // new rows
  unsigned nCols = inMat.cols(); // new cols
  gLinear::gRowMatrix<AssignType> outMat(nRows, nCols);
  for (unsigned r = 0; r < nRows; r++) { // new indices
    for (unsigned c = 0; c < nCols; c++) {
      outMat(r,c) = inMat(nRows-1-r, nCols-1-c);
    } // for
  } // for
  return outMat;
} // fliplrud

} // namespace prlite
