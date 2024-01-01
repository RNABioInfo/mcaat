//test myLib.h
#include <iostream>
#include <string>
#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "cycle_finder.h"

using namespace std;

int main(int argc, char** argv) {
    //use output_handler
    std::string sdbg_file = "/vol/d/proof_of_concept/data/" + std::string(argv[1]) + "/graph/graph";
    int length_bound = atoi(argv[2]);

    SDBG sdbg;
    sdbg.LoadFromFile(sdbg_file.c_str());
    
    cout << "Cycle Algorithm Start" << endl;

    // Fix: Find out why not output doesn't work
    CycleFinder cycle_finder(sdbg, length_bound, 47, std::string(argv[1]).c_str());
    cout << "Cycle Algorithm End" << endl;
    return 0;
}
