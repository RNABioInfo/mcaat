//test myLib.h
#include <iostream>
#include <string>
#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "cycle_finder.h"
#include "helpers/system_resources.h"

using namespace std;

int main(int argc, char** argv) {
    //use output_handler
    std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[1]) + "/graph/graph";
    int length_bound = atoi(argv[2]);

    SDBG sdbg;
    sdbg.LoadFromFile(sdbg_file.c_str());
    
    cout << "Cycle Algorithm Start" << endl;
    SystemResources system_resources;
    system_resources.AskUser(); 
    
    CycleFinder cycle_finder(sdbg, length_bound, 47, std::string(argv[1]).c_str());
    cout << "Cycle Algorithm End" << endl;
    return 0;
}
