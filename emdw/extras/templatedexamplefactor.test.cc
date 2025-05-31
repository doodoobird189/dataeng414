/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Unit tests for TemplatedExampleFactor

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


//#ifndef NO_TESTING

// patrec headers
#include "error.hpp"
#include "testing.hpp"  // TestCase

// emdw headers
#include "emdw.hpp"
#include "templatedexamplefactor.hpp"


// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <limits>

using namespace std;
using namespace emdw;

class BaseTestTemplatedExampleFactor : public TestCase
{
public:

  /** setup to be run before each test. */
  virtual void setup()
    {

    }


  /** clean up after. */
  virtual void tearDown()
    {

    }

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

};

class TestTemplatedExampleFactor : public BaseTestTemplatedExampleFactor
{
public:

  /** a name. */
  const string& isA()
    {
      static const string name = "TestTemplatedExampleFactor";
      return name;
    }


  /** run the tests. */
  virtual void runTests() {
    LOGDEBUG1(pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

  try {

    // #################### Add your testing code here ####################

    // #####################################
    // ### define some generic variables ###
    // #####################################

    TemplatedExampleFactor<unsigned> tef;

    // #####################
    // ### CLASS MEMBERS ###
    // #####################


    // ********************
    // *** Constructors ***
    // ********************

    // *****************************
    // *** Polymorphic functions ***
    // *****************************

    // *********************
    // *** Other members ***
    // *********************


    // ########################
    // ### FACTOR OPERATORS ###
    // ########################

    // **************
    // *** sample ***
    // **************

    // *****************
    // *** normalize ***
    // *****************


    // **************
    // *** absorb ***
    // **************


    // **************
    // *** cancel ***
    // **************


    // *******************
    // *** marginalize ***
    // *******************


    // **********************
    // *** observe/reduce ***
    // **********************

    // *********************
    // *** inplaceDampen ***
    // *********************


    // #################### END OF TESTS ##################################

  } // try

  catch (char msg[]) {
    LOGERROR(pLog, "Testing", "TemplatedExampleFactor threw an exception!");
    THROW(TestFailed, msg);
  } // catch

  // catch (char const* msg) {
  //   LOGERROR(pLog, "Testing", "TemplatedExampleFactor threw an exception!");
  //   THROW(TestFailed, msg);
  // } // catch

  catch (const string& msg) {
    LOGERROR(pLog, "Testing", "TemplatedExampleFactor threw an exception!");
    THROW(TestFailed, msg);
  } // catch

  catch (const exception& e) {
    LOGERROR( pLog, "Testing", "TemplatedExampleFactor threw an exception!");
    LOGERROR( pLog, "Testing", e.what() );
    throw(e);
  } // catch

  catch(...) {
    LOGERROR(pLog, "Testing", "TemplatedExampleFactor threw an exception!");
    THROW(TestFailed, "");
  } // catch


  } // runTests
};

REGISTER_TESTCASE(TestTemplatedExampleFactor);
