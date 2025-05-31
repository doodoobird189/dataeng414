/*
  $Id$
*/

#ifndef PRLITE_TESTING_HPP
#define PRLITE_TESTING_HPP

#include "prlite_logging.hpp"
#include "prlite_my_any.hpp"

#include <list>
#include <map>
#include <cmath>  // trunc, log10
#include <limits>  // numeric_limits

namespace prlite{

extern std::map<std::string, Any> testConfig;

class TestFailed : public std::exception {
  std::string msg;
public:
  TestFailed(const std::string& _m = "TestFailed: unknown cause.") throw()
    : msg(_m) { }
  TestFailed(const TestFailed& e) throw()
    : std::exception(e), msg(e.msg) { }
  virtual ~TestFailed() throw() { }
  virtual const char* what() const throw() { return msg.c_str(); }
};

/** This is a simple testsuite.  Look at myclass.test.cc for an
    example/template to use. */
class TestCase {
private:
  /// not copyable.
  TestCase(const TestCase&);

public:

  TestCase();

  /** Custom name to use for printing prettily. */
  virtual const std::string& isA() = 0;

  /** setup to be run before each test. */
  virtual void setup() = 0;

  /** clean up after. */
  virtual void tearDown() = 0;

  /** run the TestCase's tests. */
  virtual void runTests() = 0;


  /** run all the tests. */
  static void runAllTests();

  /**
   * add a TestCase to the internals by reference. You don't need to
   * delete it, and it must be created with new.
   */
  static void registerTestCase(TestCase* suite);

  /** Check true */
  void assertTrue(bool cond) {
    PRLITE_CHECK(cond, prlite::TestFailed, "Expected true.");
  }

  /** check false */
  void assertFalse(bool cond) {
    PRLITE_CHECK(!cond, prlite::TestFailed, "Expected false.");
  }

  /** Check equality. */
#define PRLITE_assertEqual(A,B) \
  PRLITE_CHECK(A == B, prlite::TestFailed, "Not equal: " << A << " != " << B)

  /** Check floating point equality. numbers are multiplied by mult,
      and then trunc()'d, and then checked for equality.

      NOTE: This returns true for nans!
  */
  bool almostEqual(double a, double b, unsigned decimalPlaces) {
    // check for nan
    if (a == a && b == b) {
      // not nan. check inf.
      static const double inf = std::numeric_limits<double>::infinity();
      static const double minf = -std::numeric_limits<double>::infinity();
      if ( (a == inf || a == minf) || (b == inf || b == minf) ) {
  return a == b;
      }
      else {
  // not nan or inf

  // make sure a is the smallest
  if (fabs(a)>fabs(b)) {
    return fabs(a-b) <= (fabs(b)/std::pow(10.0,static_cast<double>(decimalPlaces)));
  }
  else {
    return fabs(a-b) <= (fabs(a)/std::pow(10.0,static_cast<double>(decimalPlaces)));
  }
      }
    }
    else {
      return (!(a==a)) && (!(b==b));
    }
    // dummy for dumb compilers' warnings
    return false;
  }

  void assertAlmostEqual(double a, double b, unsigned decimalPlaces) {
    PRLITE_CHECK(almostEqual(a,b,decimalPlaces), prlite::TestFailed,
    "Not equal to within " << decimalPlaces << " decimal places: "
    << a << " != " << b);
  }


  /** Check inequality. */
  template<class A, class B>
  void assertNotEqual(const A& a, const B& b) {
    PRLITE_CHECK(a != b, prlite::TestFailed, "Equal: " << a << " == " << b);
  }

  /** Check exceptions */
#define PRLITE_assertThrows(EXC, DO) \
  try { \
    { DO } \
    PRLITE_THROW(prlite::TestFailed, "Expected exception " << #EXC  << " was not thrown."); \
  } \
  catch (const EXC&) {  }

  virtual ~TestCase();

private:
  static std::list<TestCase*>& allTests();
};  // TestCase

#define PRLITE_REGISTER_TESTCASE(TYPE)                                 \
  static TYPE tester_##TYPE

} // namespace prlite

#endif // PRLITE_TESTING_HPP
