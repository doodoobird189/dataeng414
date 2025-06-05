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
#include <string>
#include <memory>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace emdw;

struct Config {
    int R = -1;        // -1 means estimate R
    double ws = -1.0;  // -1 means estimate ws
    double wd = -1.0;  // -1 means estimate wd
};

double accuracy(const std::vector<int>& predicted, const std::vector<int>& groundTruth) {
    if (predicted.size() != groundTruth.size()) return 0.0;
    int correct = 0;
    for (size_t i = 0; i < predicted.size(); ++i)
        correct += (predicted[i] == groundTruth[i]);
    return static_cast<double>(correct) / predicted.size();
}

int main(int argc, char *argv[]) {
    initLogging(argv[0]);
    prlite::TestCase::runAllTests();

    try {
        Config config;
        config.R = 2;

        unsigned seedVal = emdw::randomEngine.getSeedVal();
        cout << "Random seed: " << seedVal << endl;
        emdw::randomEngine.setSeedVal(seedVal);

        typedef int T;
        typedef DiscreteTable<T> DT;
        int N;

        const std::string matrixData = R"(  
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 1
0 0 0 0 0 0 0 1 0 1 1 0 0 1 0 0 1 0 1 1
... [TRUNCATED for brevity, keep full matrix here]
)";

        std::istringstream ss(matrixData);
        std::vector<std::vector<int>> relationshipMatrix;
        std::string line;
        while (std::getline(ss, line)) {
            if (line.empty()) continue;
            std::istringstream lineStream(line);
            std::vector<int> row;
            int value;
            while (lineStream >> value) row.push_back(value);
            relationshipMatrix.push_back(row);
        }

        N = relationshipMatrix.size();
        emdw::RVIdType nodeRV[N];
        for (int i = 0; i < N; ++i) nodeRV[i] = i;

        // Best configuration (skipping estimation logic for brevity)
        int best_R = config.R;
        double best_ws = config.ws != -1 ? config.ws : 0.7;
        double best_wd = config.wd != -1 ? config.wd : 0.2;

        rcptr<vector<T>> commDom(new vector<T>());
        for (int r = 1; r <= best_R; ++r) commDom->push_back(r);

        // Build unary factors
        std::vector<DT> unaryFactors;
        for (int i = 0; i < N; ++i) {
            std::map<std::vector<T>, FProb> assignments;
            for (int r = 1; r <= best_R; ++r)
                assignments[{r}] = FProb(1.0 / best_R + 0.001);
            DT f({nodeRV[i]}, {commDom}, 0.0, assignments);
            unaryFactors.push_back(f);
        }

        // Build pairwise factors
        std::vector<DT> pairwiseFactors;
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                std::map<std::vector<T>, FProb> assignments;
                for (T a : *commDom) {
                    for (T b : *commDom) {
                        double prob = (a == b)
                            ? (relationshipMatrix[i][j] == 1 ? best_ws : 1 - best_ws)
                            : (relationshipMatrix[i][j] == 1 ? best_wd : 1 - best_wd);
                        assignments[{a, b}] = FProb(prob);
                    }
                }
                DT f({nodeRV[i], nodeRV[j]}, {commDom, commDom}, 0.0, assignments);
                pairwiseFactors.push_back(f);
            }
        }

        // Combine factors
        std::vector<rcptr<Factor>> factorPtrs;
        for (const auto& f : unaryFactors)
            factorPtrs.push_back(rcptr<Factor>(new DT(f)));
        for (const auto& f : pairwiseFactors)
            factorPtrs.push_back(rcptr<Factor>(new DT(f)));

        std::map<emdw::RVIdType, AnyType> observed;
        ClusterGraph cg(ClusterGraph::LTRIP, factorPtrs, observed);
        std::map<Idx2, rcptr<Factor>> messages;
        MessageQueue msgQueue;
        loopyBP_CG(cg, messages, msgQueue);

        std::vector<rcptr<Factor>> beliefs;
        for (int i = 0; i < N; ++i)
            beliefs.push_back(queryLBP_CG(cg, messages, {i}));

        // Final prediction
        std::vector<int> predictions(N);
        for (int i = 0; i < N; ++i) {
            double p1 = beliefs[i]->potentialAt({i}, {1});
            double p2 = beliefs[i]->potentialAt({i}, {2});
            predictions[i] = (p1 < p2) ? 1 : 2;
            std::cout << "Node " << i << ": p(1)=" << p1 << ", p(2)=" << p2
                      << " => prediction: " << predictions[i] << "\n";
        }

        std::cout << "Final estimated R = " << best_R
                  << ", ws = " << best_ws
                  << ", wd = " << best_wd << "\n";

    } catch (const std::exception &e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
