char TESTING_CC_VERSION[] = "$Id$";

/**
   Author: A.S. Visagie
   Date: 2006/05/11
   A simple testing framework.
*/

#include "prlite_testing.hpp"
#include "prlite_error.hpp"

#include <iostream>

using namespace std;

namespace prlite{

map<string, Any> testConfig;

TestCase::TestCase()
{
  // Automatically register objects
  TestCase::registerTestCase(this);
}

TestCase::~TestCase()
{

}


void TestCase::runAllTests()
{
  list<TestCase*>::iterator iter = TestCase::allTests().begin();

  int good(0);
  bool isSetup = false;

  for ( ; iter != allTests().end(); iter++) {
    TestCase& test = **iter;
    try {
      #ifdef DEBUG_MODE
      cerr << "running " << test.isA() << ", " << flush;
      #endif

      isSetup = false;
      #ifdef DEBUG_MODE
      cerr << "setup, " << flush;
      #endif
      test.setup();
      isSetup = true;

      #ifdef DEBUG_MODE
      cerr << "test, " << flush;
      #endif
      test.runTests();

      #ifdef DEBUG_MODE
      cerr << "tear-down, " << flush;
      #endif
      test.tearDown();

      #ifdef DEBUG_MODE
      cerr << "OK" << endl;
      #endif

      ++good;
    }
    catch (PatRecError& e) {
      cerr << endl << "  " << "Caught PatRecError: " << e.what() << endl;
    }
    catch (TestFailed& e) {
      cerr << endl << "  " << e.what() << endl;
    }
    catch (exception& e) {
      cerr << endl << "  " << "Caught exception: " << e.what() << endl;
    }

    if (isSetup) {
      test.tearDown();
    }
  }

  cerr << endl << endl
       << "+------------------------" << endl
       << "| passed: " << good << endl
       << "| failed: " << TestCase::allTests().size() - good << endl
       << "+------------------------" << endl;
}


void TestCase::registerTestCase(TestCase* suite)
{
  allTests().push_back(suite);
}


list<TestCase*>& TestCase::allTests()
{
  static list<TestCase*> __allTests;
  return __allTests;
}

} // namespace prlite
