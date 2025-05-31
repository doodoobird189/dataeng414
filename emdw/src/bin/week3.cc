/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_logging.hpp" // initLogging
#include "prlite_testing.hpp"

// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"

// standard headers
#include <iostream> // cout, endl, flush, cin, cerr
#include <cctype>   // toupper
#include <string>   // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <random>

using namespace std;
using namespace emdw;

// ##################################################################
//  Some example code. To compile this, go to the emdw/build
//  directory and do a:
//  cmake ../; make -j7 example
//  To run this while in the build directory, do a:
//  src/pmr/example
//
//  For your own stuff, make a copy of this one to start with. Then
//  edit the CMakeLists.txt (also in this directory) by adding your
//  new target in the same way as this example.
// ##################################################################

int main(int, char *argv[])
{

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try
  {

    //*********************************************************
    // Some random generator seeding. Just keep this as is
    //*********************************************************

    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout << seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    //*********************************************************
    // Predefine some types and constants
    //*********************************************************

    typedef int T;               // The type of the values that the RVs can take on
    typedef DiscreteTable<T> DT; // DT now is a short-hand for DiscreteTable<int>
    double defProb = 0.0;        // Any unspecified probs will default to this.
    rcptr<vector<T>> montyDom(   // Lists the values that a particular RV can take on
        new vector<T>{0, 1});

    //*********************************************************
    // Define the RVs
    //*********************************************************

    // The enum statement here predefines two RV ids: the id of X is 0
    // and the id of Y is 1. This is easy enough in very simple
    // problems, for more complex situations involving many RVs this
    // becomes cumbersome and we will need a datastructure such as a
    // map to save the RV ids in. Consult the userguide for more on
    // this.

    enum
    {
      b0,
      b1,
      b2,
      b3,
      b4,
      b5,
      b6,
      r0,
      r1,
      r2,
      r3,
      r4,
      r5,
      r6
    };

    //*********************************************************
    // Set up a discrete factor (in several ways) over two binary
    // RVs specifying that they must have odd parity (i.e. their
    // values will always differ).
    //***************************************************

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // The most direct declaration. We show this as an example of
    // construction with a basic set of parameters. See the class
    // specific constructor from line 109 in
    // src/emdw-factors/discretetable.hpp for more detail on the
    // exact types of each variable.
    //
    // IMPORTANT: However, you will instead use a dynamic
    // declaration (lower down), because that will allow you to
    // access via its abstract category namely a Factor.
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    rcptr<Factor> b4given = uniqptr<DT>(new DT(
        {b0, b1, b2, b4},                         // The ids of the two variables (ascending
                                                  // order). If not presorted, the order of the
                                                  // variables will get re-arranged in the class.
        {montyDom, montyDom, montyDom, montyDom}, // The domains over which they can vary
        defProb,                                  // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0, 0, 0}, (1)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 0, 1, 1}, (1)}, // and the 0.5 is its probabilty.
            {{0, 1, 0, 1}, (1)},
            {{0, 1, 1, 0}, (1)},
            {{1, 0, 0, 1}, (1)},
            {{1, 0, 1, 0}, (1)},
            {{1, 1, 0, 0}, (1)},
            {{1, 1, 1, 1}, (1)},
        }));

    rcptr<Factor> b5given = uniqptr<DT>(new DT(
        {b0, b2, b3, b5},                         // The ids of the two variables (ascending
                                                  // order). If not presorted, the order of the
                                                  // variables will get re-arranged in the class.
        {montyDom, montyDom, montyDom, montyDom}, // The domains over which they can vary
        defProb,                                  // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0, 0, 0}, (1)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 0, 1, 1}, (1)}, // and the 0.5 is its probabilty.
            {{0, 1, 0, 1}, (1)},
            {{0, 1, 1, 0}, (1)},
            {{1, 0, 0, 1}, (1)},
            {{1, 0, 1, 0}, (1)},
            {{1, 1, 0, 0}, (1)},
            {{1, 1, 1, 1}, (1)},
        }));

    rcptr<Factor> b6given = uniqptr<DT>(new DT(
        {b0, b1, b3, b6},                         // The ids of the two variables (ascending
                                                  // order). If not presorted, the order of the
                                                  // variables will get re-arranged in the class.
        {montyDom, montyDom, montyDom, montyDom}, // The domains over which they can vary
        defProb,                                  // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0, 0, 0}, (1)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 0, 1, 1}, (1)}, // and the 0.5 is its probabilty.
            {{0, 1, 0, 1}, (1)},
            {{0, 1, 1, 0}, (1)},
            {{1, 0, 0, 1}, (1)},
            {{1, 0, 1, 0}, (1)},
            {{1, 1, 0, 0}, (1)},
            {{1, 1, 1, 1}, (1)},
        }));

    std::cout << __FILE__ << __LINE__ << ": " << *b4given << std::endl; // displays the factor
    std::cout << __FILE__ << __LINE__ << ": " << *b5given << std::endl; // displays the factor
    std::cout << __FILE__ << __LINE__ << ": " << *b6given << std::endl
              << std::endl
              << std::endl; // displays the factor

    rcptr<Factor> result = b4given->absorb(b5given)->absorb(b6given);
    // rcptr<Factor> finalResult = result->absorb(b6given);
    result = result->normalize();

    std::cout << __FILE__ << __LINE__ << ": " << *result << std::endl
              << std::endl
              << std::endl;

    std::cout << __FILE__ << __LINE__ << ": " << *(result->observeAndReduce({b0, b1, b2, b3}, {1, 0, 1, 0})->normalize()) << std::endl;

    // R0 - R6 STUFF

    // question 3.3 e)

    rcptr<Factor> r0given = uniqptr<DT>(new DT(
        {b0, r0},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r1given = uniqptr<DT>(new DT(
        {b1, r1},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r2given = uniqptr<DT>(new DT(
        {b2, r2},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r3given = uniqptr<DT>(new DT(
        {b3, r3},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r4given = uniqptr<DT>(new DT(
        {b4, r4},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r5given = uniqptr<DT>(new DT(
        {b5, r5},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

    rcptr<Factor> r6given = uniqptr<DT>(new DT(
        {b6, r6},             // The ids of the two variables (ascending
                              // order). If not presorted, the order of the
                              // variables will get re-arranged in the class.
        {montyDom, montyDom}, // The domains over which they can vary
        defProb,              // The default probability for unspecified allocations
        {
            // The explicitly specified probabilities:
            {{0, 0}, (0.9)}, // Note: the {0,1} is the allocation to {X,Y},
            {{0, 1}, (0.1)}, // and the 0.5 is its probabilty.
            {{1, 0}, (0.1)},
            {{1, 1}, (0.9)},
        }));

        // question f
        std::cout << "QUESTION 3.3F: " << std::endl;
        std::cout << *(r0given->observeAndReduce({b0, r0}, {1,1}));
        std::cout << *(r1given->observeAndReduce({b1, r1}, {1,1}));
        std::cout << *(r2given->observeAndReduce({b2, r2}, {1,1}));
        std::cout << *(r3given->observeAndReduce({b3, r3}, {1,1}));
        std::cout << *(r4given->observeAndReduce({b4, r4}, {1,1}));
        std::cout << *(r5given->observeAndReduce({b5, r5}, {1,1}));
        std::cout << *(r6given->observeAndReduce({b6, r6}, {1,0}));

        // question 3.3 g

        rcptr<Factor> complete_result = result->absorb(r0given)->absorb(r1given)->absorb(r2given)->absorb(r3given)->absorb(r4given)->
        absorb(r5given)->absorb(r6given);
        complete_result = complete_result->normalize();

        std::cout << "PRINTING JOING DISTRIBUTION OF P(B0...R0):" << std::endl;
        std::cout << *complete_result << std::endl << std::endl << endl;

        std::cout << "PRINTING POSTERIOR MARGINAL BELIEFS FOR EACH OF THE BITS" << std::endl << std::endl;

        std::cout << "PRINTING FOR B0=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b0}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B1=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b1}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B2=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b2}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B3=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b3}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B4=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b4}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B5=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b5}, {1})->normalize()) << std::endl << std::endl;

        std::cout << "PRINTING FOR B6=1:" << std::endl;
        std::cout << __FILE__ << __LINE__ << ": " << *(complete_result->observeAndReduce({b6}, {1})->normalize()) << std::endl << std::endl;

// from emdw folder run this: cd build; cmake ../; make -j7 week3; src/bin/week3; cd ..
        

    return 0; // tell the world that all is fine
  } // try

  catch (char msg[])
  {
    cerr << msg << endl;
  } // catch

  // catch (char const* msg) {
  //   cerr << msg << endl;
  // } // catch

  catch (const string &msg)
  {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const exception &e)
  {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  } // catch

  catch (...)
  {
    cerr << "An unknown exception / error occurred\n";
    throw;
  } // catch

} // main
