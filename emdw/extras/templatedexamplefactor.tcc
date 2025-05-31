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
          factor. Then replace "TemplatedExampleFactor" with your intended new
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
#include "sortindices.hpp"
#include "vecset.hpp"

// standard headers
#include <new>  // placement new
#include <math.h>
#include <random>

using namespace std;
using namespace emdw;

/**************************The TemplatedExampleFactor class members***********************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper
// function. Boilerplate.
template <typename AssignType>
void TemplatedExampleFactor<AssignType>::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new TemplatedExampleFactor_Marginalize<AssignType>() );
  } // if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new TemplatedExampleFactor_InplaceNormalize<AssignType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new TemplatedExampleFactor_Normalize<AssignType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new TemplatedExampleFactor_InplaceAbsorb<AssignType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new TemplatedExampleFactor_Absorb<AssignType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new TemplatedExampleFactor_InplaceCancel<AssignType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new TemplatedExampleFactor_Cancel<AssignType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new TemplatedExampleFactor_ObserveAndReduce<AssignType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper = uniqptr<InplaceDamper>( new TemplatedExampleFactor_InplaceDamping<AssignType>() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new TemplatedExampleFactor_Sampler<AssignType>() );
  } // if

} // supplementWithDefaultOperators

// default constructor
template <typename AssignType>
TemplatedExampleFactor<AssignType>::TemplatedExampleFactor(
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
template <typename AssignType>
TemplatedExampleFactor<AssignType>::TemplatedExampleFactor(const TemplatedExampleFactor<AssignType>& p)
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

template <typename AssignType>
TemplatedExampleFactor<AssignType>::TemplatedExampleFactor(
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

} // class specific constructor

template <typename AssignType>
TemplatedExampleFactor<AssignType>::~TemplatedExampleFactor() {
  // TIP: Typically empty, except if your class contains raw pointers
} // destructor

//========================= Operators and Conversions =========================

template <typename AssignType>
bool TemplatedExampleFactor<AssignType>::operator==(const TemplatedExampleFactor<AssignType>& p) const {
  if (vars != p.vars) {return false;} // if

  // YCH

} // operator==

// Boilerplate
template <typename AssignType>
bool TemplatedExampleFactor<AssignType>::operator!=(const TemplatedExampleFactor<AssignType>& p) const {
  return !operator==(p);
} // operator !=

// Boilerplate code
template <typename AssignType>
bool TemplatedExampleFactor<AssignType>::isEqual(const Factor* rhsPtr) const {
  const TemplatedExampleFactor<AssignType>* dwnPtr = dynamic_cast<const TemplatedExampleFactor<AssignType>*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// Boilerplate
template <typename AssignType>
unsigned TemplatedExampleFactor<AssignType>::classSpecificConfigure(
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

  this->~TemplatedExampleFactor<AssignType>(); // Destroy existing

  new(this) TemplatedExampleFactor<AssignType>(
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

template <typename AssignType>
std::ostream& TemplatedExampleFactor<AssignType>::txtWrite(std::ostream& file) const {

  // YCH Remember the factor operators

  return file;
} // txtWrite

template <typename AssignType>
std::istream& TemplatedExampleFactor<AssignType>::txtRead(std::istream& file) {

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
template <typename AssignType>
emdw::RVVals TemplatedExampleFactor<AssignType>::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
template <typename AssignType>
uniqptr<Factor> TemplatedExampleFactor<AssignType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
template <typename AssignType>
void TemplatedExampleFactor<AssignType>::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> TemplatedExampleFactor<AssignType>::absorb(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
template <typename AssignType>
void TemplatedExampleFactor<AssignType>::inplaceAbsorb(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> TemplatedExampleFactor<AssignType>::cancel(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
template <typename AssignType>
void TemplatedExampleFactor<AssignType>::inplaceCancel(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> TemplatedExampleFactor<AssignType>::marginalize(const emdw::RVIds& variablesToKeep,bool presorted,
                                                       FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> TemplatedExampleFactor<AssignType>::observeAndReduce(const emdw::RVIds& variables,
                                                            const emdw::RVVals& assignedVals,
                                                            bool presorted,
                                                            FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
template <typename AssignType>
double TemplatedExampleFactor<AssignType>::inplaceDampen(const Factor* rhsPtr,
                                                double df,
                                                FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

template <typename AssignType>
TemplatedExampleFactor<AssignType>* TemplatedExampleFactor<AssignType>::copy(
  const emdw::RVIds& newVars,
  bool presorted) const {

  if ( newVars.size() ) { // replace the variables

    // YCH

  } // if newVars

  // identical copy
  return new TemplatedExampleFactor<AssignType>(*this);

} // copy

template <typename AssignType>
TemplatedExampleFactor<AssignType>*
TemplatedExampleFactor<AssignType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  // YCH

} // vacuousCopy

template <typename AssignType>
unsigned TemplatedExampleFactor<AssignType>::noOfVars() const {
  return vars.size();
} // noOfVars

template <typename AssignType>
emdw::RVIds TemplatedExampleFactor<AssignType>::getVars() const {
  return vars;
} // getVars

template <typename AssignType>
emdw::RVIdType TemplatedExampleFactor<AssignType>::getVar(unsigned varNo) const {
  return vars[varNo];
} // getVar

//====================== Other Virtual Member Functions =======================

template <typename AssignType>
double TemplatedExampleFactor<AssignType>::potentialAt(
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

template <typename AssignType>
double TemplatedExampleFactor<AssignType>::distance(const Factor* rhsPtr) const {

  const TemplatedExampleFactor<AssignType>* dwnPtr = dynamic_cast<const TemplatedExampleFactor<AssignType>*>(rhsPtr);
  ASSERT(dwnPtr, "TemplatedExampleFactor<AssignType>::distance : rhs is not a TemplatedExampleFactor");
  ASSERT(vars == dwnPtr->vars, "The variables do not match");

  // YCH

} // distance

//======================= Parameterizations (non-virtual) =====================


//======================= Non-virtual Member Functions ========================

//================================== Friends ==================================

/// output. Boilerplate
template<typename T>
std::ostream& operator<<( std::ostream& file,
                          const TemplatedExampleFactor<T>& n ){
  return n.txtWrite(file);
} // operator<<

/// input. Boilerplate
template<typename T>
std::istream& operator>>( std::istream& file,
                          TemplatedExampleFactor<T>& n ){
  return n.txtRead(file);
} // operator>>


//======================= Factor operators ====================================


// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample   **************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& TemplatedExampleFactor_Sampler<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_Sampler<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
emdw::RVVals TemplatedExampleFactor_Sampler<AssignType>::sample(const TemplatedExampleFactor<AssignType>* lhsPtr) {

  // YCH

} // sample


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& TemplatedExampleFactor_InplaceNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_InplaceNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
void TemplatedExampleFactor_InplaceNormalize<AssignType>::inplaceProcess(TemplatedExampleFactor<AssignType>* lhsPtr) {

  TemplatedExampleFactor<AssignType>& lhs(*lhsPtr);

  // YCH

} // inplaceProcess

// Identify object

template <typename AssignType>
const std::string& TemplatedExampleFactor_Normalize<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_Normalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boilerplate
template <typename AssignType>
Factor* TemplatedExampleFactor_Normalize<AssignType>::process(const TemplatedExampleFactor<AssignType>* lhsPtr) {

  // make a copy of the factor
  TemplatedExampleFactor<AssignType>* fPtr = new TemplatedExampleFactor<AssignType>(*lhsPtr);

  // inplace normalize it
  TemplatedExampleFactor_InplaceNormalize<AssignType> ipNorm;
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
template <typename AssignType>
const std::string& TemplatedExampleFactor_InplaceAbsorb<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_InplaceAbsorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * This operator inplace absorbs another TemplatedExampleFactor message, ie
 * does an inplace factor product.
 */
template <typename AssignType>
void TemplatedExampleFactor_InplaceAbsorb<AssignType>::inplaceProcess(TemplatedExampleFactor<AssignType>* lhsPtr,
                                      const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate TemplatedExampleFactor<AssignType>*
  const TemplatedExampleFactor<AssignType>* rhsPtr = dynamic_cast<const TemplatedExampleFactor<AssignType>*>(rhsFPtr);
  ASSERT(rhsPtr, " TemplatedExampleFactor_InplaceAbsorb::inplaceProcess : rhs not a valid TemplatedExampleFactor : " << typeid(*rhsFPtr).name() );
  const TemplatedExampleFactor<AssignType>& rhs(*rhsPtr);
  TemplatedExampleFactor<AssignType>& lhs(*lhsPtr);

  // varU - the union of vars and rhs.vars
  // u2l  - the corresponding indices in union that came from lhs.vars
  // u2r  - the corresponding indices in union that came from rhs.vars
  std::vector<int> u2l, u2r;
  emdw::RVIds varU = sortedUnion(lhs.vars, rhs.vars, u2l, u2r);
  size_t newDim = varU.size();

  // YCH

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& TemplatedExampleFactor_Absorb<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_Absorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boilerplate
template <typename AssignType>
Factor* TemplatedExampleFactor_Absorb<AssignType>::process(const TemplatedExampleFactor<AssignType>* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  TemplatedExampleFactor<AssignType>* fPtr = new TemplatedExampleFactor<AssignType>(*lhsPtr);

  // inplace normalize it
  TemplatedExampleFactor_InplaceAbsorb<AssignType> ipAbsorb;
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
template <typename AssignType>
const std::string& TemplatedExampleFactor_InplaceCancel<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_InplaceCancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * This operator inplace cancels another TemplatedExampleFactor message, ie does an
 * inplace factor division. It is assumed that we will only cancel
 * messages that were previously absorbed. If not our result might go
 * negative definite.
 */
template <typename AssignType>
void TemplatedExampleFactor_InplaceCancel<AssignType>::inplaceProcess(TemplatedExampleFactor<AssignType>* lhsPtr,
                                      const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate TemplatedExampleFactor<AssignType>*
  const TemplatedExampleFactor<AssignType>* rhsPtr = dynamic_cast<const TemplatedExampleFactor<AssignType>*>(rhsFPtr);
  ASSERT(rhsPtr, " TemplatedExampleFactor_InplaceAbsorb::inplaceProcess : rhs not a valid TemplatedExampleFactor : " << typeid(*rhsFPtr).name() );
  const TemplatedExampleFactor<AssignType>& rhs(*rhsPtr);
  TemplatedExampleFactor<AssignType>& lhs(*lhsPtr);

  std::vector<int> l2i; // which indices in lhs are shared
  std::vector<int> r2i; // with corresponding indices in rhs
  emdw::RVIds varI = sortedIntersection(lhs.vars, rhs.vars, l2i, r2i);  //intersection of vars
  ASSERT(rhs.vars.size() == varI.size(), "Factor division only supports full sharing of rhs variables");
  size_t newDim = lhs.vars.size();

  // YCH

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& TemplatedExampleFactor_Cancel<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_Cancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boilerplate
template <typename AssignType>
Factor* TemplatedExampleFactor_Cancel<AssignType>::process(const TemplatedExampleFactor<AssignType>* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  TemplatedExampleFactor<AssignType>* fPtr = new TemplatedExampleFactor<AssignType>(*lhsPtr);

  // inplace normalize it
  TemplatedExampleFactor_InplaceCancel<AssignType> ipCancel;
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
template <typename AssignType>
const std::string& TemplatedExampleFactor_Marginalize<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_Marginalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* TemplatedExampleFactor_Marginalize<AssignType>::process(const TemplatedExampleFactor<AssignType>* lhsPtr,
                                            const emdw::RVIds& toKeep,
                                            bool presorted) {

  const TemplatedExampleFactor<AssignType>& lhs(*lhsPtr);
  int origDim = lhs.vars.size();
  ASSERT(toKeep.size(), "TemplatedExampleFactor_Marginalize : Can not marginalize everything out");

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
  ASSERT(theVars.size(), "TemplatedExampleFactor_Marginalize : Can not marginalize everything out");
  ASSERT(!rni.size(), "TemplatedExampleFactor_Marginalize : Can not marginalize to variables not in scope");

  // YCH

} // process


/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& TemplatedExampleFactor_ObserveAndReduce<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_ObserveAndReduce<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* TemplatedExampleFactor_ObserveAndReduce<AssignType>::process(const TemplatedExampleFactor<AssignType>* lhsPtr,
                                          const emdw::RVIds& theVars,
                                          const emdw::RVVals& anyVals,
                                          bool presorted) {


  // and a convenient reference for the lhs factor
  const TemplatedExampleFactor<AssignType>& lhs(*lhsPtr);
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
template <typename AssignType>
const std::string& TemplatedExampleFactor_InplaceDamping<AssignType>::isA() const {
  static const std::string CLASSNAME("TemplatedExampleFactor_InplaceDamping<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and old
// parameters - this linear combination is determined by the df parameter
template <typename AssignType>
double
TemplatedExampleFactor_InplaceDamping<AssignType>::inplaceProcess(
  TemplatedExampleFactor<AssignType>* newPtr,
  const Factor* oldPtr,
  double oldWeight) { // oldWeight is the weight for rhs

  // A typedef for the actual internal assignment vector type
  //typedef std::vector<AssignType> std::vector<AssignType>;

  TemplatedExampleFactor<AssignType>& lhs(*newPtr);

  // YCH

} // inplaceProcess
