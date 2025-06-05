/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"

// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"
#include "textblockio.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <random>
#include <sstream>
#include <vector>
#include <fstream>  // ifstream, ofstream

using namespace std;
using namespace emdw;

//##################################################################
// Some example code. To compile this, go to the emdw/build
// directory and do a:
// cmake ../; make -j7 example
// To run this while in the build directory, do a:
// src/pmr/example
//
// For your own stuff, make a copy of this one to start with. Then
// edit the CMakeLists.txt (also in this directory) by adding your
// new target in the same way as this example.
//##################################################################

// Function to count lines in a file (number of nodes)
int countLines(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }
    return count(istreambuf_iterator<char>(file), 
                istreambuf_iterator<char>(), '\n');
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
        int value;
        int community = -1;
        
        // Find which column has the 1 (community assignment)
        int col = 0;
        while (iss >> value) {
            if (value == 1) {
                community = col + 1; // Convert to 1-based index
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

int main(int, char *argv[]) {

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {

    //*********************************************************
    // Some random generator seeding. Just keep this as is
    //*********************************************************

    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout <<  seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    gLinear::gRowMatrix<int> adjecency_matrix;
    adjecency_matrix = loadBlock<int>("testing.txt");
    std::cout << "Adjacency matrix loaded with " 
              << adjecency_matrix.rows() << " nodes and "
              << adjecency_matrix.cols() << " features." << std::endl;

    return 0; // tell the world that all is fine
  } // try

  catch (char msg[]) {
    cerr << msg << endl;
  } // catch

  // catch (char const* msg) {
  //   cerr << msg << endl;
  // } // catch

  catch (const string& msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const exception& e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  } // catch

  catch(...) {
    cerr << "An unknown exception / error occurred\n";
    throw;
  } // catch

} // main
