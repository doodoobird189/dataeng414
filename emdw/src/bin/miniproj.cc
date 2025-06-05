// patrec headers
#include "prlite_logging.hpp" // initLogging
#include "prlite_testing.hpp"

// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"

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
#include <sstream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
using namespace emdw;

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

    int R = 2; // number of communities
    int N;  // number of nodes
    
    // Setting values for ws and wd
    double ws = 0.7;
    double wd = 0.2;

    double esti_ws = 0.9;
    double esti_wd = 0.1;

    typedef int T;                             // The type of the values that the RVs can take on
    typedef DiscreteTable<T> DT;               // DT now is a short-hand for DiscreteTable<int>
    double defProb = 0.0;                      // Any unspecified probs will default to this.
    rcptr<vector<T>> commDom(new vector<T>()); // Lists the values that a particular RV can take on

    for (int i = 1; i <= R; i++)
    {
      commDom->push_back(i);
    }

    // enter input data - make sure there are no unnecessary spaces
    const std::string matrixData = R"(  
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 1
0 0 0 0 0 0 0 1 0 1 1 0 0 1 0 0 1 0 1 1
0 0 0 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1
0 0 1 0 1 0 1 1 1 0 0 1 1 0 1 1 0 0 0 1
0 0 1 1 0 1 1 0 1 0 1 1 1 0 1 1 1 0 0 0
0 0 1 0 1 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0
0 0 1 1 1 1 0 1 1 1 1 0 1 0 1 1 0 1 1 1
1 1 1 1 0 0 1 0 1 1 0 0 1 1 1 0 1 0 1 0
0 0 1 1 1 0 1 1 0 1 0 1 0 0 1 1 0 0 1 1
0 1 1 0 0 0 1 1 1 0 0 0 0 1 0 0 1 0 0 0
0 1 1 0 1 0 1 0 0 0 0 1 1 0 1 1 0 0 1 1
0 0 1 1 1 0 0 0 1 0 1 0 1 0 1 0 0 0 1 1
0 0 1 1 1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1
0 1 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0
0 0 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 0 1
0 0 1 1 1 1 1 0 1 0 1 0 1 0 1 0 0 0 1 1
1 1 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0
1 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0
0 1 1 0 0 1 1 1 1 0 1 1 1 0 0 1 0 0 0 0
1 1 1 1 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 0
)";

    // Populating the relationship matrix with the values from the input matrix.
    // The relationship matrix is the edge matrix
    std::istringstream ss(matrixData);
    std::vector<std::vector<int>> relationshipMatrix;
    std::string line;

    while (std::getline(ss, line)) {
        std::istringstream lineStream(line);
        std::vector<int> row;
        int value;
        while (lineStream >> value) {
            row.push_back(value);
        }

        if (!row.empty()) {
            relationshipMatrix.push_back(row);
        }
    }

    // determining number of nodes in given data i.e N
    std::istringstream iss(matrixData);
    line = "";
    int n = 0;

    // Skip empty lines and find first valid line
    while (std::getline(iss, line)) {
        if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) {
            std::istringstream linestream(line);
            std::string val;
            while (linestream >> val) ++n;
            break;
        }
    }
    N = n;

    // Creating random variables for the nodes
    emdw::RVIdType nodeRV[N];
    for (int i = 0; i < N; i++)
    {
      nodeRV[i] = i;
    }

    std::srand(42);

    // zi on its own
    // Create a factor for each node with uniform probability 1.0/R + e on each state
    std::vector<DT> factorList1;
    for (int i = 0; i < N; i++)
    {
      std::map<std::vector<T>, FProb> assignments;
      for (int val = 1; val <= R; val++)
      {
        // FProb is assumed to have .prob and an optional flag (here set to zero)
        assignments[{val}] = FProb((1.0 / R) + 0.0001 + (static_cast<double>(rand()) / RAND_MAX) * (0.001 - 0.0001));  // noise is added to the uniform distribution to break symmetry
      }
      DT newFactor({nodeRV[i]}, {commDom}, 0.0, assignments);
      factorList1.push_back(newFactor);
    }

    // zj on its own
    // create a factor for each node with uniform probability 1.0/R on each state
    std::vector<DT> factorList2;
    for (int i = 0; i < N; i++)
    {
      std::map<std::vector<T>, FProb> assignments;
      for (int val = 1; val <= R; val++)
      {
        // FProb is assumed to have .prob and an optional flag (here set to zero)
        assignments[{val}] = FProb(1.0 / R);
      }
      DT newFactor({nodeRV[i]}, {commDom}, 0.0, assignments);
      factorList2.push_back(newFactor);
    }

    // Create factors between every unique pair of nodes (i, j with i < j)
    std::vector<DT> pairFactors;
    for (int i = 0; i < N; i++)
    {
      for (int j = i + 1; j < N; j++)
      {
        std::map<std::vector<T>, FProb> assignments;
        // Iterate over each possible assignment to the two nodes.
        for (auto a : *commDom)
        {
          for (auto b : *commDom)
          {
            double prob;
            if (a == b)
            {
              // Same value
              prob = (relationshipMatrix[i][j] == 1) ? ws : (1 - ws);
            }
            else
            {
              // Different values
              prob = (relationshipMatrix[i][j] == 1) ? wd : (1 - wd);
            }
            assignments[{a, b}] = FProb(prob);
          }
        }
        // Create a vector of domain pointers for the two variables.
        std::vector<rcptr<vector<T>>> pairDomains = {commDom, commDom};
        DT newFactor({nodeRV[i], nodeRV[j]}, pairDomains, 0.0, assignments);
        pairFactors.push_back(newFactor);
      }
    }

    // The node factors are added to a list of pair factors
    for (auto dt: factorList1)
    {
      pairFactors.push_back(dt);
    }

    // converting pairfactors to rcptr<factor> so that we can use them in the cluster graph
    std::vector<rcptr<Factor>> factorPtrs;
    for (const auto &dt : pairFactors)
    {
      factorPtrs.push_back(rcptr<Factor>(new DT(dt)));
    }

    // You may optionally supply an observation map (empty if none).
    std::map<emdw::RVIdType, AnyType> observed; // empty observations in this example.

    ClusterGraph cg(ClusterGraph::LTRIP, factorPtrs, observed);//, true);

    // Optionally, export the graph to visualize it
    cg.exportToGraphViz("myClusterGraph");

    std::cout << "ClusterGraph created and exported to myClusterGraph.dot" << std::endl;

    /*
        std::cout << pairFactors[1] << std::endl;
        std::cout << relationshipMatrix[0][1] << std::endl;
        std::cout << "I love cake" << std::endl;
        // Assuming that RVIds is a typedef for a container of ints (or similar)
    auto rvids = pairFactors[1].getVars();
    std::cout << "Random Variables in factor 1: ";
    for (auto rv : rvids) {
        std::cout << rv << " ";
    }
        */

    std::cout << factorList1[0] << std::endl;
    std::cout << pairFactors[0] << std::endl;
    //std::cout << pairFactors[0].get << std::endl;

    /* LOOPY BELIEF PROPAGATION */

    // ... [code that creates your factors and builds a ClusterGraph 'cg' as already shown] ...

    // Set up the observation map (empty if you have no evidence)

    // Create containers for messages and the update queue.
    // (Idx2 is used to index messages between clusters; MessageQueue is your scheduler type.)
    std::map<Idx2, rcptr<Factor>> messages;
    MessageQueue msgQueue;

    // Set LBP parameters
    double dampingFactor = 0.2;    // controls how much new messages affect the current messages
    double deltaThresh = 1e-4;     // change threshold to decide convergence of message updates
    unsigned maxNoOfAbsorbs = 100; // maximum absorb operations allowed (will be scaled internally)

    // Run loopy belief propagation on the cluster graph.
    unsigned absorbCount = loopyBP_CG(cg,
                                      messages,
                                      msgQueue);
    std::cout << "Loopy BP converged after " << absorbCount << " absorb operations." << std::endl;

    // Now query a belief (i.e. marginal) on your cluster graph.
    // For example, let’s query the marginal for a particular variable or a set of variables.
    // Assume that emdw::RVIds is a typedef (e.g., to std::vector<int>) of random variable IDs.
    emdw::RVIds queryVars = {3};
    rcptr<Factor> belief = queryLBP_CG(cg, messages, queryVars);

    // Assuming your Factor class or derived DiscreteTable supports printing,
    // you might then output the result.
    std::cout << "Belief for variable(s) ";
    for (auto v : queryVars)
    {
      std::cout << v << " ";
    }
    std::cout << ":\n";
    belief->txtWrite(std::cout);
    std::cout << std::endl;

    // EM algorithm to estimate the parameters: ws and wd
    vector<rcptr<Factor>> node_beliefs;
    for (int i = 0; i<N; i++) 
    {
      node_beliefs.push_back(queryLBP_CG(cg, messages, {i}));
    }

  
    
    double ws_numer = 0.0;
    double ws_denom = 0.0;
    double wd_numer = 0.0;
    double wd_denom = 0.0;
    double p_same = 0.0;
    double p_diff = 0.0;
    for (int i = 0; i < N; i++)
    {
      for (int j = i + 1; j < N; j++) // sum over numerator and denominator
      {
        p_same = 0.0;
        p_diff = 0.0;

        // Calculate the probabilities for the same and different states
        // using the node beliefs and the relationship matrix.
        // The potentialAt method is assumed to return the probability of a given state.
        for (int r = 1; r <= R; r++) // sum over p same
        {
          p_same += node_beliefs[i]->potentialAt({i}, {r}) * node_beliefs[j]->potentialAt({j}, {r});
        }
        p_diff = 1 - p_same;

        ws_numer += p_same * relationshipMatrix[i][j];
        ws_denom += p_same;

        wd_numer += p_diff * relationshipMatrix[i][j];
        wd_denom += p_diff;
      }
    }
    esti_ws = ws_numer / ws_denom;
    esti_wd = wd_numer / wd_denom;

    std::cout << "Estimated ws: " << esti_ws << std::endl;
    std::cout << "Estimated wd: " << esti_wd << std::endl




    /* Loopy Belief Update */

    // --- Setup for Loopy BP ---
    // Create containers for messages and the message update queue.
    // LBP parameters: damping factor, convergence threshold, and maximum operations.
    dampingFactor = 0.2;  // How strongly new messages affect the current ones.
    deltaThresh = 1e-4;   // Convergence threshold for message updates.
    maxNoOfAbsorbs = 100; // Upper bound on absorb operations (scaled internally).

    std::map<Idx2, rcptr<Factor>> messages2;
    MessageQueue msgQueue2;
    // Run loopy belief update on the cluster graph.
    absorbCount = loopyBU_CG(cg,
                             messages2,
                             msgQueue2);
    std::cout << "Loopy BU converged after " << absorbCount << " absorb operations." << std::endl;

    // --- Querying the Beliefs ---
    // Define which variables you want to obtain the marginal for.
    // For instance, if you want the marginal for variable 19:
    //queryVars = {2};

    // Obtain the belief (marginal probability factor) for the specified query variables.
    belief = queryLBU_CG(cg, messages, queryVars);

    // Output the belief. Here, we assume that the Factor (or DiscreteTable)
    // has a txtWrite method that outputs to an ostream.
    std::cout << "Belief for variable(s) ";
    for (auto v : queryVars)
      std::cout << v << " ";
    std::cout << ":\n";
    belief->txtWrite(std::cout);
    std::cout << std::endl;

    //std::cout << pairFactors[2] << std::endl;
    //std::cout << "Size: " << pairFactors.size() << std::endl;
    //std::cout << belief->getProb({1}) << std::endl;
    std::cout << belief->potentialAt({3},{2}) << std::endl;

    std::cout << std::endl;
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
