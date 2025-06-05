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

/********************************************************************************************************
 * CONFIGURATION BELOW
 * You can set the values of R, ws, and wd here. If you set them to -1, they will be estimated.
 * You can also set the file names for the adjacency matrix and ground truth communities.
 ********************************************************************************************************
*/

struct Config {
  int R = -1;         // -1 means estimate R
  double ws = -1;   // -1 means estimate ws
  double wd = -1;   // -1 means estimate wd
};

// change the names of the files below to match your data files
const string adjecencyMatrixFile = "A5.txt";
const string groundTruthFile = "V5_ground_truth.txt";

std::string trim(const std::string& str) {
    size_t start = 0;
    while (start < str.size() && std::isspace(str[start])) ++start;
    size_t end = str.size();
    while (end > start && std::isspace(str[end - 1])) --end;
    return str.substr(start, end - start);
}

double accuracy(const std::vector<int>& predicted, const std::vector<int>& groundTruth) {
    if (predicted.size() != groundTruth.size()) {
        std::cerr << "Error: Vectors are of different sizes." << std::endl;
        return 0.0;
    }

    int correct = 0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == groundTruth[i]) {
            ++correct;
        }
    }

    return static_cast<double>(correct) / predicted.size();
}

// Function to read adjacency matrix with dynamic size
vector<vector<int>> readAdjacencyMatrix(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening adjacency matrix file: " << filename << endl;
        exit(1);
    }

    vector<vector<int>> adjMatrix;
    string line;
    
    while (getline(file, line)) {
        vector<int> row;
        istringstream iss(line);
        int value;
        while (iss >> value) {
            row.push_back(value);
        }
        adjMatrix.push_back(row);
    }

    // Validate square matrix
    size_t numNodes = adjMatrix.size();
    for (const auto& row : adjMatrix) {
        if (row.size() != numNodes) {
            cerr << "Error: Adjacency matrix is not square!" << endl;
            exit(1);
        }
    }

    return adjMatrix;
}

// Function to read ground truth communities
vector<int> readGroundTruth(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening ground truth file: " << filename << endl;
        exit(1);
    }

    vector<int> communities;
    string line;

    while (getline(file, line)) {
        istringstream iss(line);
        double value;
        int community = -1;
        int col = 0;

        while (iss >> value) {
            if (fabs(value - 1.0) < 1e-6) {  // tolerate floating point precision
                community = col + 1; // 1-based indexing
                break;
            }
            col++;
        }

        if (community == -1) {
            cerr << "Error: No community assignment found for node" << endl;
            exit(1);
        }

        communities.push_back(community);
    }

    return communities;
}

void copyGroundTruthFile(const std::string& sourcePath) {
    std::ifstream source(sourcePath, std::ios::binary);
    std::ofstream dest("community_truth.txt", std::ios::binary);
    
    if (!source.is_open()) {
        std::cerr << "Error: Could not open source file " << sourcePath << std::endl;
        return;
    }
    if (!dest.is_open()) {
        std::cerr << "Error: Could not create community_truth.txt" << std::endl;
        return;
    }

    dest << source.rdbuf(); // Copy file contents

    /*
    // Verify copy was successful
    if (source && dest) {
        std::cout << "Successfully copied ground truth file" << std::endl;
    } else {
        std::cerr << "Error occurred during file copy" << std::endl;
    } */
}


int main(int argc, char *argv[])
{
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {
    Config config;

    // Set manually here if desired
    // uncomment the following lines to set fixed values ---- leave commented to estimate
    //config.R = 3;     // Comment this line to estimate R
    //config.ws = 0.7;  // Comment this to estimate ws
    //config.wd = 0.2;  // Comment this to estimate wd

    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout << seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    typedef int T;
    typedef DiscreteTable<T> DT;
    int N;

    // READING DATA
    std::vector<std::vector<int>> relationshipMatrix = readAdjacencyMatrix(adjecencyMatrixFile);
    std::vector<int> community_assignments = readGroundTruth(groundTruthFile);
    N = relationshipMatrix.size();


/*
    const std::string matrixData = R"(  
0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0
1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0
1 1 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0
1 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 1
0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1
0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 1 0 1
0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1 0
)";

// enter gournd truths data - make sure there are no unnecessary spaces
    const std::string groundTruthData = R"(
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
1 0
0 1
1 0
1 0
1 0
1 0
0 1
0 1
1 0
1 0
0 1
1 0
0 1
1 0
0 1
0 1
0 1
0 1
0 1
0 1
0 1
0 1
0 1
0 1
0 1
0 1
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
*/
    emdw::RVIdType nodeRV[N];
    for (int i = 0; i < N; i++) nodeRV[i] = i;
    std::srand(42);

    int best_R = config.R;
    double best_ws = config.ws;
    double best_wd = config.wd;
    vector<rcptr<Factor>> best_node_beliefs;

    // Determine what needs to be estimated
    bool estimate_R = (config.R == -1);
    bool estimate_ws = (config.ws == -1);
    bool estimate_wd = (config.wd == -1);
    if (estimate_ws)
    {
        config.ws = 0.8;
    }
    if(estimate_wd)
    {
        config.wd = 0.2;
    }

    if (estimate_R || estimate_ws || estimate_wd) {
      struct ModelResult {
        int R;
        double log_likelihood;
        double ws, wd;
        vector<rcptr<Factor>> beliefs;
      };

      const int R_min = estimate_R ? 2 : config.R;
      const int R_max = estimate_R ? 6 : config.R;
      std::vector<ModelResult> modelResults;

      for (int R_try = R_min; R_try <= R_max; ++R_try) {
        rcptr<vector<T>> commDom(new vector<T>());
        for (int r = 1; r <= R_try; ++r) commDom->push_back(r);

        std::vector<DT> factorList1;
        for (int i = 0; i < N; ++i) {
          std::map<std::vector<T>, FProb> assignments;
          for (int val = 1; val <= R_try; val++) {
            assignments[{val}] = FProb((1.0 / R_try) + 0.001 + (static_cast<double>(rand()) / RAND_MAX) * 0.01);
          }
          DT newFactor({nodeRV[i]}, {commDom}, 0.0, assignments);
          factorList1.push_back(newFactor);
        }

        double esti_ws = estimate_ws ? 0.9 : config.ws;
        double esti_wd = estimate_wd ? 0.1 : config.wd;
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
          if (!estimate_ws && !estimate_wd) break; // skip updates if not estimating either
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

        modelResults.push_back({R_try, log_likelihood, esti_ws, esti_wd, node_beliefs});
      }

      auto best = *max_element(modelResults.begin(), modelResults.end(), [](const ModelResult& a, const ModelResult& b) {
        return a.log_likelihood < b.log_likelihood;
      });

      best_R = best.R;
      best_ws = best.ws;
      best_wd = best.wd;
      best_node_beliefs = best.beliefs;
    }

    // If no estimation was done, or beliefs weren't computed:
    if (best_node_beliefs.empty()) {
      rcptr<vector<T>> commDom(new vector<T>());
      for (int r = 1; r <= best_R; ++r) commDom->push_back(r);

      std::vector<DT> factorList1;
      for (int i = 0; i < N; ++i) {
        std::map<std::vector<T>, FProb> assignments;
        for (int val = 1; val <= best_R; val++) {
          assignments[{val}] = FProb((1.0 / best_R) + 0.001 + (static_cast<double>(rand()) / RAND_MAX) * 0.01);
        }
        DT newFactor({nodeRV[i]}, {commDom}, 0.0, assignments);
        factorList1.push_back(newFactor);
      }

      std::vector<DT> pairFactors;
      for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
          std::map<std::vector<T>, FProb> assignments;
          for (auto a : *commDom) {
            for (auto b : *commDom) {
              double prob = (a == b)
                ? ((relationshipMatrix[i][j] == 1) ? best_ws : (1 - best_ws))
                : ((relationshipMatrix[i][j] == 1) ? best_wd : (1 - best_wd));
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

      best_node_beliefs.clear();
      for (int i = 0; i < N; ++i)
        best_node_beliefs.push_back(queryLBP_CG(cg, messages, {i}));
    }

    std::cout << "\nFinal node beliefs:\n";
    for (int i = 0; i < N; ++i) {
      std::cout << "Node " << i << ": ";
      best_node_beliefs[i]->txtWrite(std::cout);
      std::cout << std::endl;
    }
    std::cout << "\nFinal estimated R = " << best_R << ", ws = " << best_ws << ", wd = " << best_wd << "\n";

    /*
    for (int i = 0; i < N; i++) {
        std::vector<double> probs;
        for(int r = 1; r<best_R+1; r++)
        {
            probs.push_back(best_node_beliefs[i]->potentialAt({i}, {r}));
        }
        auto minIt = std::min_element(probs.begin(), probs.end());
        int minIndex = std::distance(probs.begin(), minIt);
        predictions[i] = minIndex + 1; // Adjust for 0-based index
    }
        
    for (int i = 0; i < N; i++)
    {

        if (best_node_beliefs[i]->potentialAt({i}, {1}) < best_node_beliefs[i]->potentialAt({i}, {2}))
        {
            predictions[i] = 1;
            std::cout << best_node_beliefs[i]->potentialAt({i}, {1}) << " > " << best_node_beliefs[i]->potentialAt({i}, {2}) << std::endl;
        }
        if (best_node_beliefs[i]->potentialAt({i}, {1}) > best_node_beliefs[i]->potentialAt({i}, {2}))
        {
            predictions[i] = 2;
            std::cout << best_node_beliefs[i]->potentialAt({i}, {1}) << " <= " << best_node_beliefs[i]->potentialAt({i}, {2}) << std::endl;
        }
    }

    */
   /*
    for (int i = 0; i < N; i++)
    {
        double max_potential = -std::numeric_limits<double>::infinity();
        int best_label = -1;

        // Assuming labels are 0-based and go up to R-1
        for (int label = 0; label < best_R; label++)
        {
            double potential = best_node_beliefs[i]->potentialAt({i}, {label});
            if (potential > max_potential)
            {
                max_potential = potential;
                best_label = label;
            }
        }

        predictions[i] = best_label;
        std::cout << "Node " << i << " best label: " << best_label
                  << " with potential: " << max_potential << std::endl;
    }
*/
    vector<int> predictions(N);
    for(int n=0; n<N; n++)
    {
         vector<double> probs;
        for(int r=1; r<best_R+1; r++)
        {
            probs.push_back(best_node_beliefs[n]->potentialAt({n}, {r}));
        }
        double max = probs[0];
        for(int i=1; i<probs.size(); i++)
        {
            if (probs[i] > max)
            {
                max = probs[i];
            }
        }
        for(int i=0; i<probs.size(); i++)   // find the max
        {
            if (probs[i] == max)
            {
                predictions[n] = i+1;
            }
        }
        //std::cout << "Node " << n << " best label: " << predictions[n]
                  //<< " with potential: " << max << std::endl;
    }
/*
    // reading ground truth data for accuracy calculation
    // Accuracy of the predictions
    std::vector<int> community_assignments(N);  // Fixed size allocation
    std::istringstream stream(groundTruthData);
    std::string line2;
    int nodeIndex = 0;

    while (std::getline(stream, line2) && nodeIndex < N) {
    line2 = trim(line2);
    if (line2.empty()) {
        continue;  // skip empty lines but don't increment nodeIndex
    }

    std::istringstream lineStream(line2);
    std::vector<double> probabilities;
    double value;
    while (lineStream >> value) {
        probabilities.push_back(value);
    }

    if (!probabilities.empty()) {
        auto maxIt = std::max_element(probabilities.begin(), probabilities.end());
        int community = static_cast<int>(std::distance(probabilities.begin(), maxIt)) + 1;
        community_assignments[nodeIndex] = community;
        nodeIndex++;  // Only increment here for valid lines
    }
}
*/
    // QUANTITATIVE OUTPUT - FINAL PREDICTIONS
    std::cout << "Final predictions for each node: \n";
    std::cout << "Format is as follows: Node <node_index>: <predicted_community> ---- <actual_community>\n";
    for(int i = 0; i < N; i++)
    {
        std::cout << "Node " << i << ": " << predictions[i] << " ---- " << community_assignments[i] << std::endl;
    }
    std::cout << std::endl;
    
    /*
    // Output the community assignments
    std::cout << "Community Assignments:\n";
    for (int c : community_assignments) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
*/
    std::cout << "Accuracy: " << accuracy(predictions, community_assignments) * 100 << "%\n";

  // QUALITATIVE OUTPUT
  // save predictions to a file and later use it for visualisation in a python script
std::ofstream outfile("community_predictions.txt");

// Write node index and predicted community
for(int i = 0; i < N; i++) {
    outfile << i << " " << predictions[i] << "\n";
}
outfile.close();

// Also save the adjacency matrix for graph structure
std::ofstream adjfile("adjacency_matrix.txt");
for(const auto& row : relationshipMatrix) {
    for(int val : row) adjfile << val << " ";
    adjfile << "\n";
}
adjfile.close();

// Copy ground truth file for reference
copyGroundTruthFile(groundTruthFile);

    return 0;
  } catch (const exception &e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  } catch (...) {
    cerr << "An unknown exception / error occurred\n";
    throw;
  }
}
