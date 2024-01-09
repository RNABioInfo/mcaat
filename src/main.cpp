//test myLib.h
#include <iostream>
#include <string>
#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "cycle_finder.h"
#include "helpers/system_resources.h"
#include "graph_writer.h"

using namespace std;



int main(int argc, char** argv) {
    //use output_handler
        //use output_handler
    if (argc == 4 && std::string(argv[1]) == "u_mode"){
        std::cout << "Usage: ./cycle_finder <mode> <genome_name> <length_bound>" << std::endl;
        std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[2]) + "/graph/graph";
        int length_bound = atoi(argv[3]);

        SDBG sdbg;
        sdbg.LoadFromFile(sdbg_file.c_str());
        //GraphWriter gw(sdbg);
        //gw.WriteNodes();
        cout << "Cycle Algorithm Start" << endl;
        
        CycleFinder cycle_finder(sdbg, length_bound, 47, std::string(argv[2]).c_str());
        cout << "Cycle Algorithm End" << endl;
    }
    else{
        std::cout << "Usage: ./cycle_finder <mode> <genome_name> <length_bound>" << std::endl;
    }
    return 0;
}
