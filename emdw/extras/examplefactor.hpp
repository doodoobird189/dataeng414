#ifndef EXAMPLEFACTOR_HPP
#define EXAMPLEFACTOR_HPP

/*******************************************************************************

          AUTHORS: JA du Preez

          DATE: September 2018

          PURPOSE: A code template to facilitate implementing a new
          Factor type.

          Start by copying this {.hpp,.cc} file to a name that is the
          full lower case version of the name for your new
          factor. Then replace "ExampleFactor" with your intended new
          factor name. As per our convention for classes, your new
          factor name should start with an uppercase letter. Various
          code fragments that probably will prove handy are retained
          in the following - at least consider them carefully before
          changing it. In some places code will be marked as
          "boiler-plate" your probably should not tamper with that.

          You also need to implement a examplefactor.test.cc file
          containing the unit tests for your factor and
          operators. This is not supplied here. Use one of the
          existing .test.cc's as a template.

          To be able to compile, update the src/CMakeLists.txt file to
          include this new code. Do this and check that you can
          compile.

          Now the actual work starts, supplying the implementation
          detail for your new classes. Watch out for:

           TIP -> A comment to take note of
           YCH -> Your Code Here

          It is important to first familiarise yourself with the
          syntax and semantics of the various factor and factor
          operator members described in factor.hpp and
          factoroperator.hpp. You HAVE to abide by their rules.

          Please provide decent descriptions for the interfaces of your
          functions, it is time-consuming to first study the code
          before one understands the parameters of a function.

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "emdw.hpp"
#include "factor.hpp" // Factor
#include "factoroperator.hpp" // FactorOperator
#include "anytype.hpp"

// patrec headers

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

// Some forward declarations
class ExampleFactor_InplaceNormalize;
class ExampleFactor_InplaceMaxNormalize;
class ExampleFactor_InplaceAbsorb;
class ExampleFactor_Absorb;
class ExampleFactor_InplaceCancel;
class ExampleFactor_Marginalize;
class ExampleFactor_MaxMarginalize;
class ExampleFactor_ObserveAndReduce;
class ExampleFactor_InplaceDamping;
class ExampleFactor_Sampler;

/******************************************************************************/
/************************ ExampleFactor **************************************/
/******************************************************************************/

/**
 * TIP: Give a decent description of the class and its purpose here.
 */
class ExampleFactor : public Factor {

  // TIP: these are some friend declarations, basically because it is
  // expedient to allow factor operators direct access to the private
  // data of its associated factor.
  friend class ExampleFactor_InplaceNormalize;
  friend class ExampleFactor_InplaceMaxNormalize;
  friend class ExampleFactor_InplaceAbsorb;
  friend class ExampleFactor_Absorb;
  friend class ExampleFactor_InplaceCancel;
  friend class ExampleFactor_Marginalize;
  friend class ExampleFactor_MaxMarginalize;
  friend class ExampleFactor_ObserveAndReduce;
  friend class ExampleFactor_InplaceDamping;
  friend class ExampleFactor_Sampler;

  //============================ Traits and Typedefs ============================
private:

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:

  /**
   * @brief Default constructor.
   *
   * @param theVars Each variable in the PGM will be identified
   * with a specific integer that uniquely ids it.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  ExampleFactor(
    const emdw::RVIds& theVars = {},
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  ///
  ExampleFactor(const ExampleFactor& p);

  /// TIP: If you have pointers as data members in your class, you
  /// probably also need a move construcor. Below we assume this is
  /// not the case.
  ExampleFactor(ExampleFactor&& st) = default;

  /**
    * @brief Class Specific Constructor TIP: This is the
    * detailed/primary constructor for the class. Make sure you
    * describe its parameters thoroughly, not only the syntax but also
    * semantics/ intent. The tmpVal parameter is bogus, we only
    * include it to differentiate from the default constructor.
    *
    */
  ExampleFactor(
    const emdw::RVIds& theVars,
    unsigned tmpVal,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  ///
  virtual ~ExampleFactor();

  //========================= Operators and Conversions =========================
public:

  /// TIP: If you have pointers of some type in the class, you should
  /// supply these operators. Below we assume this to not be the case.
  ExampleFactor& operator=(const ExampleFactor& d) = default;
  ExampleFactor& operator=(ExampleFactor&& d) = default;

  ///
  bool operator==(const ExampleFactor& d) const;

  ///
  bool operator!=(const ExampleFactor& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /// TIP: this one is a run-time equivalent of the class specific
  /// constructor. It basically rebuilds an existing object.
  unsigned classSpecificConfigure(
    const emdw::RVIds& theVars,
    unsigned tmpval,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  /// input from stream
  virtual std::istream& txtRead(std::istream& file);
  /// output to stream
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from Factor (at the top-level). They return uniqptr so
  /// that we can chain them. Boilerplate.

  //------------------Family 0: eg Sample

  virtual emdw::RVVals sample(Sampler* procPtr) const;

  //------------------Family 1: eg Normalize

  virtual void inplaceNormalize(FactorOperator* procPtr = 0);
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr = 0) const;

  //------------------Family 2: eg Absorb,Cancel

  virtual void inplaceAbsorb(const Factor* rhsPtr,
                            FactorOperator* procPtr = 0);
  /// result will inherit all the operators from the lhs factor
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  virtual void inplaceCancel(const Factor* rhsPtr,
                            FactorOperator* procPtr = 0);
  /// result will inherit all the operators from the lhs factor
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  //------------------Family 3 eg Marginalize

  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                     bool presorted = false,
                                     FactorOperator* procPtr = 0) const;

  //------------------Family 4: eg Observe, ObserveAndReduce

  virtual uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                           const emdw::RVVals& assignedVals,
                                           bool presorted = false,
                                           FactorOperator* procPtr = 0) const;

  //------------------Family 5: eg inplaceDampen

  /**
   * NOTE the parameter is the OLD message and will be weighted with df
   */
  virtual double inplaceDampen(const Factor* oldMsg,
                               double df,
                               FactorOperator* procPtr = 0);

  //===================== Required Virtual Member Functions =====================
public:

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace, one-for-one AND IN
   * THIS ORDER, the original SORTED set.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual ExampleFactor* copy(
    const emdw::RVIds& newVars = {},
    bool presorted = false) const;

  /**
   * make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   */
  virtual ExampleFactor* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /// note rhs is a raw pointer
  bool isEqual(const Factor* rhsPtr) const;

  /// TIP: This default implementation should do the trick, but do
  /// provide your own it if necessary.
  double distanceFromVacuous() const {return Factor::distanceFromVacuous();}

  ///
  virtual unsigned noOfVars() const;

  /// returns the ids of the variables the factor is defined on
  virtual emdw::RVIds getVars() const;

  ///
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //====================== Other Virtual Member Functions =======================
public:

  /// returns the factor potential at a given variable assignment
  virtual double potentialAt(const emdw::RVIds& theVars,
                             const emdw::RVVals& theirVals,
                             bool presorted = false) const;

  /**
   * Measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference. Typically based
   * on KL divergence.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Parameterizations (non-virtual) =====================
public:

  //======================= Non-virtual Member Functions ========================
public:

  //================================== Friends ==================================

public:

  /// output
  friend std::ostream& operator<<( std::ostream& file,
                                   const ExampleFactor& n );

  /// input
  friend std::istream& operator>>( std::istream& file,
                                   ExampleFactor& n );

  //=============================== Data Members ================================
private:

  emdw::RVIds vars;

  // TIP/YCH: the factor private data goes here

  // The default set of FactorOperators in use. Boilerplate.
  rcptr<Marginalizer> marginalizer;
  rcptr<InplaceNormalizer> inplaceNormalizer;
  rcptr<Normalizer> normalizer;
  rcptr<InplaceAbsorber> inplaceAbsorber;
  rcptr<Absorber> absorber;
  rcptr<InplaceCanceller> inplaceCanceller;
  rcptr<Canceller> canceller;
  rcptr<ObserveAndReducer> observeAndReducer;
  rcptr<InplaceDamper> inplaceDamper;
  rcptr<Sampler> sampler;

  //=============================== DEPRECATED ================================

}; // ExampleFactor

// TIP: Try to stick to the factor operator naming convention used
// below. It makes life easier.

// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a ExampleFactor.
 */
class ExampleFactor_Sampler: public SampleOperator< ExampleFactor > {
public:
  const std::string& isA() const;
  ///
  emdw::RVVals sample(
    const ExampleFactor* lhsPtr);
}; // ExampleFactor_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

/// This operator does an inplace normalization of a ExampleFactor.
class ExampleFactor_InplaceNormalize: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(ExampleFactor* lstPtr);
}; // ExampleFactor_InplaceNormalize

/// This operator returns a normalized version of the input factor.
class ExampleFactor_Normalize: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  Factor* process(const ExampleFactor* lstPtr);
}; // ExampleFactor_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

/// This operator does an inplace factor product with another
/// ExampleFactor.
class ExampleFactor_InplaceAbsorb: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(ExampleFactor* lhsPtr,
                      const Factor* rhsFPtr);
}; // ExampleFactor_InplaceAbsorb

/// This operator returns a factor product with another
/// ExampleFactor.
class ExampleFactor_Absorb: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  Factor* process(const ExampleFactor* lhsPtr,
                  const Factor* rhsFPtr);
}; // ExampleFactor_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

/// This operator does an inplace factor division with another
/// ExampleFactor.
class ExampleFactor_InplaceCancel: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(ExampleFactor* lhsPtr,
                      const Factor* rhsFPtr);
}; // ExampleFactor_InplaceCancel

/// This operator returns a factor division with another
/// ExampleFactor.
class ExampleFactor_Cancel: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  Factor* process(const ExampleFactor* lhsPtr,
                  const Factor* rhsFPtr);
}; // ExampleFactor_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

/** This operator returns a marginalized factor.
 */
class ExampleFactor_Marginalize: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  Factor* process(const ExampleFactor* lhsPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false);
}; // ExampleFactor_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

/**
 * After observing some variables, the scope of the resultant factor
 * reduces to ONLY that of the unobserved variables, making further
 * computations potentially much more efficient.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 */
class ExampleFactor_ObserveAndReduce: public Operator1< ExampleFactor > {
public:
  const std::string& isA() const;
  Factor* process(const ExampleFactor* lhsPtr,
                const emdw::RVIds& variables,
                const emdw::RVVals& assignedVals,
                bool presorted);
}; // ExampleFactor_ObserveAndReduce

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************

/**
 * Dampens/smooths two corresponding ExampleFactors and returns the
 * total distance between them.
 *
 * new = lhs*(1-oldWeight) + rhs*oldWeight
 */
class ExampleFactor_InplaceDamping: public Operator1< ExampleFactor > {
public:

  const std::string& isA() const;

  /**
   * @param newPtr is the new message
   *
   * @param oldPtr is the old message
   *
   * @param oldWeight is the weighting/damping_factor for the old message
   */
  double
  inplaceProcess(
    ExampleFactor* newPtr,
    const Factor* oldPtr,
    double oldWeight);

}; // ExampleFactor_InplaceDamping


#endif // EXAMPLEFACTOR_HPP
