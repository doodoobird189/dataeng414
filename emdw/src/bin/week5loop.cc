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
#include <vector>

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

    rcptr<Factor> result = b4given->absorb(b5given)->absorb(b6given);

    result = result->normalize();

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

        rcptr<Factor> b0distri = uniqptr<DT>(new DT (
            {b0},
            {montyDom},
            defProb,
            {
                {{0}, (0.5)},
                {{1}, (0.5)},
            }
        ));

        rcptr<Factor> b1distri = uniqptr<DT>(new DT (
            {b1},
            {montyDom},
            defProb,
            {
                {{0}, (0.5)},
                {{1}, (0.5)},
            }
        ));
        rcptr<Factor> b2distri = uniqptr<DT>(new DT (
            {b2},
            {montyDom},
            defProb,
            {
                {{0}, (0.5)},
                {{1}, (0.5)},
            }
        ));
        rcptr<Factor> b3distri = uniqptr<DT>(new DT (
            {b3},
            {montyDom},
            defProb,
            {
                {{0}, (0.5)},
                {{1}, (0.5)},
            }
        ));

    /*
        DEFINING CLUSTER NODES
    */
   
    //rcptr<Factor> clusternode1 = b0distri->absorb(b1distri)->absorb(b2distri)->absorb(b3distri)->normalize();
    rcptr<Factor> clusternode1 = b4given->absorb(r0given->observeAndReduce({r0}, {1}))->absorb(r1given->observeAndReduce({r1}, {0}))->
    absorb(r2given->observeAndReduce({r2}, {0}))->absorb(r4given->observeAndReduce({r4}, {0}));

    rcptr<Factor> clusternode2 = b5given->absorb(r3given->observeAndReduce({r3}, {0}))->absorb(r5given->observeAndReduce({r5}, {0}));

    rcptr<Factor> clusternode3 = b6given->absorb(r6given->observeAndReduce({r6}, {1}));


    // defining messages
    rcptr<Factor> m12 = uniqptr<DT>(new DT (
        {b0, b2, b3},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));

    rcptr<Factor> m21 = uniqptr<DT>(new DT (
        {b0, b2, b3},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));

    rcptr<Factor> m13 = uniqptr<DT>(new DT (
        {b0, b1, b3},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));

    rcptr<Factor> m31 = uniqptr<DT>(new DT (
        {b0, b1, b3},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));

    rcptr<Factor> m14 = uniqptr<DT>(new DT (
        {b0, b1, b2},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));

    rcptr<Factor> m41= uniqptr<DT>(new DT (
        {b0, b1, b2},
        {montyDom, montyDom, montyDom},
        defProb,
        {
            {{0, 0, 0}, (0.125)},
            {{0, 0, 1}, (0.125)},
            {{0, 1, 0}, (0.125)},
            {{0, 1, 1}, (0.125)},
            {{1, 0, 0}, (0.125)},
            {{1, 0, 1}, (0.125)},
            {{1, 1, 0}, (0.125)},
            {{1, 1, 1}, (0.125)},
        }
    ));


    // septets
    rcptr<Factor> sep1= uniqptr<DT>(new DT (
        {b0, b2},
        {montyDom, montyDom},
        defProb,
        {
            {{0, 0}, (0.25)},
            {{0, 1}, (0.25)},
            {{1, 0}, (0.25)},
            {{1, 1}, (0.25)}
        }
    ));

    rcptr<Factor> sep2 = uniqptr<DT>(new DT (
        {b3},
        {montyDom},
        defProb,
        {
            {{0}, (0.5)},
            {{1}, (0.5)}
        }
    ));

    rcptr<Factor> sep3= uniqptr<DT>(new DT (
        {b0, b1},
        {montyDom, montyDom},
        defProb,
        {
            {{0, 0}, (0.25)},
            {{0, 1}, (0.25)},
            {{1, 0}, (0.25)},
            {{1, 1}, (0.25)}
        }
    ));
/*
    // anti-clockwise
    rcptr<Factor> sep1new = clusternode1->marginalize({b0, b2});
    clusternode2 = sep1new->cancel(sep1)->absorb(clusternode2);
    sep1 = sep1new;

    rcptr<Factor> sep2new = clusternode2->marginalize({b3});
    clusternode3 = sep2new->cancel(sep2)->absorb(clusternode3);
    sep2 = sep2new;

    rcptr<Factor> sep3new = clusternode3->marginalize({b0, b1});
    clusternode1 = sep3new->cancel(sep3)->absorb(clusternode1);

    // clockwise
    sep1new = clusternode2->marginalize({b0, b2});
    clusternode1 = sep1new->cancel(sep1)->absorb(clusternode1);
    sep1 = sep1new;

    sep3new = clusternode1->marginalize({b0, b1});
    clusternode3 = sep3new->cancel(sep3)->absorb(clusternode3);
    sep3 = sep3new;

    sep2new = clusternode3->marginalize({b3});
    clusternode2 = sep2new->cancel(sep2)->absorb(clusternode2);
    sep2 = sep2new;
*/
    rcptr<Factor> sep1new, sep2new, sep3new;

    vector<int> distances = { 10000, 10000, 10000 };

    while(*std::max_element(distances.begin(), distances.end()) < 0.00000001) {
        // counter-clockwise
        sep1new = clusternode1->marginalize({b0, b2});
        clusternode2 = sep1new->cancel(sep1)->absorb(clusternode2);
        sep1 = sep1new;

        rcptr<Factor> sep2new = clusternode2->marginalize({b3});
        clusternode3 = sep2new->cancel(sep2)->absorb(clusternode3);
        sep2 = sep2new;

        rcptr<Factor> sep3new = clusternode3->marginalize({b0, b1});
        clusternode1 = sep3new->cancel(sep3)->absorb(clusternode1);

        // clockwise
        sep1new = clusternode2->marginalize({b0, b2});
        clusternode1 = sep1new->cancel(sep1)->absorb(clusternode1);
        distances[0] = sep1->distance(sep1new);
        sep1 = sep1new;

        sep3new = clusternode1->marginalize({b0, b1});
        clusternode3 = sep3new->cancel(sep3)->absorb(clusternode3);
        distances[1] = sep3->distance(sep3new);
        sep3 = sep3new;

        sep2new = clusternode3->marginalize({b3});
        clusternode2 = sep2new->cancel(sep2)->absorb(clusternode2);
        distances[2] = sep2->distance(sep2new);
        sep2 = sep2new;
    }

    //std::cout << *clusternode4->normalize() << std::endl;

    std::cout << "PRINTING MARGINAL BELIEFS:" << std::endl;
    std::cout << *clusternode1->marginalize({b0})->normalize() << std::endl;
    std::cout << *clusternode1->marginalize({b1})->normalize() << std::endl;
    std::cout << *clusternode1->marginalize({b2})->normalize() << std::endl;
    std::cout << *clusternode2->marginalize({b3})->normalize() << std::endl;
    std::cout << *clusternode1->marginalize({b4})->normalize() << std::endl;
    std::cout << *clusternode2->marginalize({b5})->normalize() << std::endl;
    std::cout << *clusternode3->marginalize({b6})->normalize() << std::endl;

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

