//test myLib.h
#include <iostream>
#include <string>
#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "cycle_finder.h"
#include "helpers/system_resources.h"

using namespace std;

int main(int argc, char** argv) {
    //use output_handler
   // std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[1]) + "/graph/graph";
    int length_bound = atoi(argv[1]);
    //cout << "SDBG file: " << sdbg_file << endl;
    vector<string> sdbg_files = {"CP071793.1"};
    for (string sdbg_file : sdbg_files) {
        SDBG sdbg;
        string file = "proof_of_concept/data/" + sdbg_file + "/graph/graph";
        sdbg.LoadFromFile(file.c_str());
        cout << "mCAAT" << endl;
        CycleFinder cycle_finder(sdbg, length_bound, 47, sdbg_file);
        cout << "mCAAT for " << sdbg_file << " finished"<< endl;
        sdbg;
    }
    return 0;
}
