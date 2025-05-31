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
#include <iostream>
#include <cctype>
#include <string>
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
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {
    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout << seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    typedef int T;
    typedef DiscreteTable<T> DT;
    int N;

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

    std::istringstream ss(matrixData);
    std::vector<std::vector<int>> relationshipMatrix;
    std::string line;
    while (std::getline(ss, line)) {
      std::istringstream lineStream(line);
      std::vector<int> row;
      int value;
      while (lineStream >> value) row.push_back(value);
      if (!row.empty()) relationshipMatrix.push_back(row);
    }

    std::istringstream iss(matrixData);
    line = "";
    int n = 0;
    while (std::getline(iss, line)) {
      if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) {
        std::istringstream linestream(line);
        std::string val;
        while (linestream >> val) ++n;
        break;
      }
    }
    N = n;

    emdw::RVIdType nodeRV[N];
    for (int i = 0; i < N; i++) nodeRV[i] = i;
    std::srand(42);

    const int R_min = 2;
    const int R_max = 6;

    struct ModelResult {
      int R;
      double log_likelihood;
      double ws, wd;
    };

    std::vector<ModelResult> modelResults;

    for (int R_try = R_min; R_try <= R_max; ++R_try) {
      std::cout << "\n>>> Trying R = " << R_try << " communities\n";

      rcptr<vector<T>> commDom(new vector<T>());
      for (int r = 1; r <= R_try; ++r) commDom->push_back(r);

      std::vector<DT> factorList1;
      for (int i = 0; i < N; ++i) {
        std::map<std::vector<T>, FProb> assignments;
        for (int val = 1; val <= R_try; val++) {
          assignments[{val}] = FProb((1.0 / R_try) + 0.0001 + (static_cast<double>(rand()) / RAND_MAX) * 0.001);
        }
        DT newFactor({nodeRV[i]}, {commDom}, 0.0, assignments);
        factorList1.push_back(newFactor);
      }

      double esti_ws = 0.9;
      double esti_wd = 0.1;
      double prev_ws = esti_ws, prev_wd = esti_wd;
      const int MAX_ITER = 20;
      const double TOL = 1e-4;
      vector<rcptr<Factor>> node_beliefs;

      for (int iter = 0; iter < MAX_ITER; ++iter) {
        std::vector<DT> pairFactors;
        for (int i = 0; i < N; ++i) {
          for (int j = i + 1; j < N; ++j) {
            std::map<std::vector<T>, FProb> assignments;
            for (auto a : *commDom) {
              for (auto b : *commDom) {
                double prob = (a == b)
                  ? ((relationshipMatrix[i][j] == 1) ? esti_ws : (1 - esti_ws))
                  : ((relationshipMatrix[i][j] == 1) ? esti_wd : (1 - esti_wd));
                assignments[{a, b}] = FProb(prob);
              }
            }
            DT newFactor({nodeRV[i], nodeRV[j]}, {commDom, commDom}, 0.0, assignments);
            pairFactors.push_back(newFactor);
          }
        }

        for (auto dt : factorList1) pairFactors.push_back(dt);
        std::vector<rcptr<Factor>> factorPtrs;
        for (const auto &dt : pairFactors)
          factorPtrs.push_back(rcptr<Factor>(new DT(dt)));

        std::map<emdw::RVIdType, AnyType> observed;
        ClusterGraph cg(ClusterGraph::LTRIP, factorPtrs, observed);
        std::map<Idx2, rcptr<Factor>> messages;
        MessageQueue msgQueue;
        loopyBP_CG(cg, messages, msgQueue);

        node_beliefs.clear();
        for (int i = 0; i < N; ++i)
          node_beliefs.push_back(queryLBP_CG(cg, messages, {i}));

        double ws_numer = 0, ws_denom = 0, wd_numer = 0, wd_denom = 0;
        for (int i = 0; i < N; ++i) {
          for (int j = i + 1; j < N; ++j) {
            double p_same = 0.0;
            for (int r = 1; r <= R_try; ++r)
              p_same += node_beliefs[i]->potentialAt({i}, {r}) * node_beliefs[j]->potentialAt({j}, {r});
            double p_diff = 1.0 - p_same;
            ws_numer += p_same * relationshipMatrix[i][j];
            ws_denom += p_same;
            wd_numer += p_diff * relationshipMatrix[i][j];
            wd_denom += p_diff;
          }
        }

        double new_ws = ws_numer / (ws_denom + 1e-10);
        double new_wd = wd_numer / (wd_denom + 1e-10);
        if (std::abs(new_ws - prev_ws) < TOL && std::abs(new_wd - prev_wd) < TOL) break;
        prev_ws = esti_ws = new_ws;
        prev_wd = esti_wd = new_wd;
      }

      double log_likelihood = 0.0;
      for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
          double p_ij = 0.0;
          for (int r = 1; r <= R_try; ++r) {
            for (int s = 1; s <= R_try; ++s) {
              double pr = node_beliefs[i]->potentialAt({i}, {r});
              double ps = node_beliefs[j]->potentialAt({j}, {s});
              double p_edge = (r == s) ? esti_ws : esti_wd;
              if (relationshipMatrix[i][j] == 1)
                p_ij += pr * ps * log(p_edge + 1e-10);
              else
                p_ij += pr * ps * log(1 - p_edge + 1e-10);
            }
          }
          log_likelihood += p_ij;
        }
      }

      modelResults.push_back({R_try, log_likelihood, esti_ws, esti_wd});
    }

    ModelResult best = *max_element(modelResults.begin(), modelResults.end(), [](const ModelResult& a, const ModelResult& b) {
      return a.log_likelihood < b.log_likelihood;
    });

    std::cout << "\nBest R = " << best.R << " with log-likelihood = " << best.log_likelihood << std::endl;
    std::cout << "Estimated ws = " << best.ws << ", wd = " << best.wd << std::endl;

    std::cout << std::endl;
    return 0;
  }

  catch (const exception &e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  }

  catch (...) {
    cerr << "An unknown exception / error occurred\n";
    throw;
  }
}
