/*******************************************************************************

          AUTHORS: JA du Preez

          DATE: September 2018

          PURPOSE: A code template to facilitate implementing a new
          Factor type.

          It is important to first familiarise yourself with the
          syntax and semantics of the various factor and factor
          operator members described in factor.hpp and
          factoroperator.hpp. You HAVE to abide by their rules.

          Start by copying this {.hpp,.cc} file to a name that is the
          full lower case version of the name for your new
          factor. Then replace "ExampleFactor" with your intended new
          factor name. As per our convention for classes, your new
          factor name should start with an uppercase letter. Various
          code fragments that probably will prove handy are retained
          in the following - at least consider them carefully before
          changing it. In some places code will be marked as
          "boiler-plate" your probably should not tamper with that.

          Please provide decent descriptions for the interfaces of your
          functions, it is time-consuming to first study the code
          before one understands the parameters of a function.

           TIP -> A comment to take note of
           YCH -> Your Code Here

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

//emdw headers
#include "emdw.hpp"
#include "examplefactor.hpp"
#include "sortindices.hpp"
#include "vecset.hpp"

// standard headers
#include <new>  // placement new
#include <math.h>
#include <random>

using namespace std;
using namespace emdw;

/**************************The ExampleFactor class members***********************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper
// function. Boilerplate.
void ExampleFactor::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new ExampleFactor_Marginalize() );
  } // if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new ExampleFactor_InplaceNormalize() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new ExampleFactor_Normalize() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new ExampleFactor_InplaceAbsorb() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new ExampleFactor_Absorb() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new ExampleFactor_InplaceCancel() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new ExampleFactor_Cancel() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new ExampleFactor_ObserveAndReduce() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper = uniqptr<InplaceDamper>( new ExampleFactor_InplaceDamping() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new ExampleFactor_Sampler() );
  } // if

} // supplementWithDefaultOperators

// default constructor
ExampleFactor::ExampleFactor(
  const emdw::RVIds& theVars,
  bool presorted,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler)
    : vars( theVars.size() ),
      marginalizer(theMarginalizer),
      inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper),
      sampler(theSampler) {

  supplementWithDefaultOperators();

  if ( theVars.size() ) {

    if (presorted or theVars.size() == 1) {
      vars = theVars;
    } // if

    else {
      std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );
      vars = extract<unsigned>(theVars, sorted);
    } // else

  } // if

  // YCH

} // default constructor

// copy constructor. We provide this to make sure we do not have gLinear deep/shallow issues
ExampleFactor::ExampleFactor(const ExampleFactor& p)
    : vars(p.vars),
      marginalizer(p.marginalizer),
      inplaceNormalizer(p.inplaceNormalizer),
      normalizer(p.normalizer),
      inplaceAbsorber(p.inplaceAbsorber),
      absorber(p.absorber),
      inplaceCanceller(p.inplaceCanceller),
      canceller(p.canceller),
      observeAndReducer(p.observeAndReducer),
      inplaceDamper(p.inplaceDamper),
      sampler(p.sampler) {
  // NOTE: We are assuming the operators are without memory, i.e. we
  // can apply the same operator to several factor instances without
  // any harm.

  // YCH

} // cpy ctor

ExampleFactor::ExampleFactor(
  const emdw::RVIds& theVars,
  unsigned tmpval,
  bool presorted,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler)
: vars(theVars.size()),
  marginalizer(theMarginalizer),
  inplaceNormalizer(theInplaceNormalizer),
  normalizer(theNormalizer),
  inplaceAbsorber(theInplaceAbsorber),
  absorber(theAbsorber),
  inplaceCanceller(theInplaceCanceller),
  canceller(theCanceller),
  observeAndReducer(theObserveAndReducer),
  inplaceDamper(theInplaceDamper),
  sampler(theSampler) {

  size_t nDims = vars.size();
  if (nDims == 1) {presorted = true;}

  supplementWithDefaultOperators();

  // YCH

} // class specific constructor Triangular CanonicalForm


ExampleFactor::~ExampleFactor() {
  // TIP: Typically empty, except if your class contains raw pointers
} // destructor

//========================= Operators and Conversions =========================

bool ExampleFactor::operator==(const ExampleFactor& p) const {
  if (vars != p.vars) {return false;} // if

  // YCH

} // operator==

// Boilerplate
bool ExampleFactor::operator!=(const ExampleFactor& p) const {
  return !operator==(p);
} // operator !=

// Boilerplate code
bool ExampleFactor::isEqual(const Factor* rhsPtr) const {
  const ExampleFactor* dwnPtr = dynamic_cast<const ExampleFactor*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// Boilerplate
unsigned ExampleFactor::classSpecificConfigure(
  const emdw::RVIds& theVars,
  unsigned tmpval,
  bool presorted,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler) {

  this->~ExampleFactor(); // Destroy existing

  new(this) ExampleFactor(
    theVars,
    tmpval,
    presorted,
    theMarginalizer,
    theInplaceNormalizer,
    theNormalizer,
    theInplaceAbsorber,
    theAbsorber,
    theInplaceCanceller,
    theCanceller,
    theObserveAndReducer,
    theInplaceDamper,
    theSampler); // and do an inplace construction

  return 1;
} // classSpecificConfigure

std::ostream& ExampleFactor::txtWrite(std::ostream& file) const {

  // YCH Remember the factor operators

  return file;
} // txtWrite

std::istream& ExampleFactor::txtRead(std::istream& file) {

  RVIds theVars; // TIP: Bogus, read this from file
  unsigned tmpval;  // TIP: Bogus, supply your own types

  // YCH Remember the factor operators

  classSpecificConfigure(theVars, tmpval);
  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
emdw::RVVals ExampleFactor::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
uniqptr<Factor> ExampleFactor::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
void ExampleFactor::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
uniqptr<Factor> ExampleFactor::absorb(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
void ExampleFactor::inplaceAbsorb(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
uniqptr<Factor> ExampleFactor::cancel(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
void ExampleFactor::inplaceCancel(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
uniqptr<Factor> ExampleFactor::marginalize(const emdw::RVIds& variablesToKeep,bool presorted,
                                                       FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
uniqptr<Factor> ExampleFactor::observeAndReduce(const emdw::RVIds& variables,
                                                            const emdw::RVVals& assignedVals,
                                                            bool presorted,
                                                            FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
double ExampleFactor::inplaceDampen(const Factor* rhsPtr,
                                                double df,
                                                FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

ExampleFactor* ExampleFactor::copy(
  const emdw::RVIds& newVars,
  bool presorted) const {

  if ( newVars.size() ) { // replace the variables

    // YCH

  } // if newVars

  // identical copy
  return new ExampleFactor(*this);

} // copy

ExampleFactor*
ExampleFactor::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  // YCH

} // vacuousCopy

unsigned ExampleFactor::noOfVars() const {
  return vars.size();
} // noOfVars

emdw::RVIds ExampleFactor::getVars() const {
  return vars;
} // getVars

emdw::RVIdType ExampleFactor::getVar(unsigned varNo) const {
  return vars[varNo];
} // getVar

//====================== Other Virtual Member Functions =======================

double ExampleFactor::potentialAt(
  const emdw::RVIds& theVars,
  const emdw::RVVals& theirVals,
  bool presorted) const{

  if (theVars.size() == 1) {presorted = true;}

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);
  ASSERT(sortedVars == vars, "The variables do not match: " << sortedVars << " vs " << vars);

    // YCH

} // potentialAt

double ExampleFactor::distance(const Factor* rhsPtr) const {

  const ExampleFactor* dwnPtr = dynamic_cast<const ExampleFactor*>(rhsPtr);
  ASSERT(dwnPtr, "ExampleFactor::distance : rhs is not a ExampleFactor");
  ASSERT(vars == dwnPtr->vars, "The variables do not match");

  // YCH

} // distance

//======================= Parameterizations (non-virtual) =====================


//======================= Non-virtual Member Functions ========================

//================================== Friends ==================================

/// output. Boilerplate
std::ostream& operator<<( std::ostream& file,
                          const ExampleFactor& n ){
  return n.txtWrite(file);
} // operator<<

/// input. Boilerplate
std::istream& operator>>( std::istream& file,
                          ExampleFactor& n ){
  return n.txtRead(file);
} // operator>>


//======================= Factor operators ====================================


// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample   **************************************
// ***************************************************************************

// Identify object
const std::string& ExampleFactor_Sampler::isA() const {
  static const std::string CLASSNAME("ExampleFactor_Sampler");
  return CLASSNAME;
} // isA

emdw::RVVals ExampleFactor_Sampler::sample(const ExampleFactor* lhsPtr) {

  // YCH

} // sample


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

// Identify object
const std::string& ExampleFactor_InplaceNormalize::isA() const {
  static const std::string CLASSNAME("ExampleFactor_InplaceNormalize");
  return CLASSNAME;
} // isA

void ExampleFactor_InplaceNormalize::inplaceProcess(ExampleFactor* lhsPtr) {

  ExampleFactor& lhs(*lhsPtr);

  // YCH

} // inplaceProcess

// Identify object

const std::string& ExampleFactor_Normalize::isA() const {
  static const std::string CLASSNAME("ExampleFactor_Normalize");
  return CLASSNAME;
} // isA

// Boilerplate
Factor* ExampleFactor_Normalize::process(const ExampleFactor* lhsPtr) {

  // make a copy of the factor
  ExampleFactor* fPtr = new ExampleFactor(*lhsPtr);

  // inplace normalize it
  ExampleFactor_InplaceNormalize ipNorm;
  try{
    ipNorm.inplaceProcess(fPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;
} // process


/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

// Identify object
const std::string& ExampleFactor_InplaceAbsorb::isA() const {
  static const std::string CLASSNAME("ExampleFactor_InplaceAbsorb");
  return CLASSNAME;
} // isA

/**
 * This operator inplace absorbs another ExampleFactor message, ie
 * does an inplace factor product.
 */
void ExampleFactor_InplaceAbsorb::inplaceProcess(ExampleFactor* lhsPtr,
                                      const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate ExampleFactor*
  const ExampleFactor* rhsPtr = dynamic_cast<const ExampleFactor*>(rhsFPtr);
  ASSERT(rhsPtr, " ExampleFactor_InplaceAbsorb::inplaceProcess : rhs not a valid ExampleFactor : " << typeid(*rhsFPtr).name() );
  const ExampleFactor& rhs(*rhsPtr);
  ExampleFactor& lhs(*lhsPtr);

  // varU - the union of vars and rhs.vars
  // u2l  - the corresponding indices in union that came from lhs.vars
  // u2r  - the corresponding indices in union that came from rhs.vars
  std::vector<int> u2l, u2r;
  emdw::RVIds varU = sortedUnion(lhs.vars, rhs.vars, u2l, u2r);
  size_t newDim = varU.size();

  // YCH

} // inplaceProcess

// Identify object
const std::string& ExampleFactor_Absorb::isA() const {
  static const std::string CLASSNAME("ExampleFactor_Absorb");
  return CLASSNAME;
} // isA

// Boilerplate
Factor* ExampleFactor_Absorb::process(const ExampleFactor* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  ExampleFactor* fPtr = new ExampleFactor(*lhsPtr);

  // inplace normalize it
  ExampleFactor_InplaceAbsorb ipAbsorb;
  try{
    ipAbsorb.inplaceProcess(fPtr,rhsFPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;

} // process


/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

// Identify object
const std::string& ExampleFactor_InplaceCancel::isA() const {
  static const std::string CLASSNAME("ExampleFactor_InplaceCancel");
  return CLASSNAME;
} // isA

/**
 * This operator inplace cancels another ExampleFactor message, ie does an
 * inplace factor division. It is assumed that we will only cancel
 * messages that were previously absorbed. If not our result might go
 * negative definite.
 */
void ExampleFactor_InplaceCancel::inplaceProcess(ExampleFactor* lhsPtr,
                                      const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate ExampleFactor*
  const ExampleFactor* rhsPtr = dynamic_cast<const ExampleFactor*>(rhsFPtr);
  ASSERT(rhsPtr, " ExampleFactor_InplaceAbsorb::inplaceProcess : rhs not a valid ExampleFactor : " << typeid(*rhsFPtr).name() );
  const ExampleFactor& rhs(*rhsPtr);
  ExampleFactor& lhs(*lhsPtr);

  std::vector<int> l2i; // which indices in lhs are shared
  std::vector<int> r2i; // with corresponding indices in rhs
  emdw::RVIds varI = sortedIntersection(lhs.vars, rhs.vars, l2i, r2i);  //intersection of vars
  ASSERT(rhs.vars.size() == varI.size(), "Factor division only supports full sharing of rhs variables");
  size_t newDim = lhs.vars.size();

  // YCH

} // inplaceProcess

// Identify object
const std::string& ExampleFactor_Cancel::isA() const {
  static const std::string CLASSNAME("ExampleFactor_Cancel");
  return CLASSNAME;
} // isA

// Boilerplate
Factor* ExampleFactor_Cancel::process(const ExampleFactor* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  ExampleFactor* fPtr = new ExampleFactor(*lhsPtr);

  // inplace normalize it
  ExampleFactor_InplaceCancel ipCancel;
  try{
    ipCancel.inplaceProcess(fPtr,rhsFPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;
} // process


/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

// Identify object
const std::string& ExampleFactor_Marginalize::isA() const {
  static const std::string CLASSNAME("ExampleFactor_Marginalize");
  return CLASSNAME;
} // isA

Factor* ExampleFactor_Marginalize::process(const ExampleFactor* lhsPtr,
                                            const emdw::RVIds& toKeep,
                                            bool presorted) {

  const ExampleFactor& lhs(*lhsPtr);
  int origDim = lhs.vars.size();
  ASSERT(toKeep.size(), "ExampleFactor_Marginalize : Can not marginalize everything out");

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &toKeep;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( toKeep, std::less<unsigned>() );
    tmpVars = extract<unsigned>(toKeep,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);

  std::vector<int> l2i; // which indices in lhs pot are shared
  std::vector<int> r2i; // with corresponding indices in rhs pot, not used here
  std::vector<int> lni; // indices in lhs that are not in intersection, i.e. to be marginalized out
  std::vector<int> rni; // should be empty, all variables should be in lhs scope
  emdw::RVIds theVars = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);
  ASSERT(theVars.size(), "ExampleFactor_Marginalize : Can not marginalize everything out");
  ASSERT(!rni.size(), "ExampleFactor_Marginalize : Can not marginalize to variables not in scope");

  // YCH

} // process


/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

// Identify object
const std::string& ExampleFactor_ObserveAndReduce::isA() const {
  static const std::string CLASSNAME("ExampleFactor_ObserveAndReduce");
  return CLASSNAME;
} // isA

Factor* ExampleFactor_ObserveAndReduce::process(const ExampleFactor* lhsPtr,
                                          const emdw::RVIds& theVars,
                                          const emdw::RVVals& anyVals,
                                          bool presorted) {


  // and a convenient reference for the lhs factor
  const ExampleFactor& lhs(*lhsPtr);
  int origDim = lhs.vars.size();
  if (theVars.size() == 1) {presorted = true;}

  // if observed nothing, return as is
  if (!theVars.size()) {
    return lhs.copy();
  } // if

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);

  // now i need to find out where the intersection between vars and
  // sortedVars lie. Any variables outside of these need not be
  // considered.
  std::vector<int> l2i; // which indices in lhs pot are observed (i.e. y)
  std::vector<int> r2i; // with corresponding indices in rhs pot
  std::vector<int> lni; // non-observed lhs variables ARE the new variable set (i.e. x)
  std::vector<int> rni; // remaining rhs variables (i.e. observed vars not part of this factor and thus ignored)
  emdw::RVIds varI = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);

  // YCH

} // process


// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************

// Identify object
const std::string& ExampleFactor_InplaceDamping::isA() const {
  static const std::string CLASSNAME("ExampleFactor_InplaceDamping");
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and old
// parameters - this linear combination is determined by the df parameter
double
ExampleFactor_InplaceDamping::inplaceProcess(
  ExampleFactor* newPtr,
  const Factor* oldPtr,
  double oldWeight) { // oldWeight is the weight for rhs

  // A typedef for the actual internal assignment vector type
  //typedef std::vector<AssignType> std::vector<AssignType>;

  ExampleFactor& lhs(*newPtr);

  // YCH

} // inplaceProcess
