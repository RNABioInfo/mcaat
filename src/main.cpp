//test myLib.h
#include <iostream>
#include "../libs/megahit/src/main_sdbg_build.cpp"
#include "core/cycle_finder.cpp"
#include "helpers/system_resources.h"
#include <fstream>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph_writer.h"
//#include "filtering_seqan.h"
//include <seqan/align.h>

///using namespace seqan;
using namespace std;
using std::istringstream;

int main(int argc, char** argv) {
    /*
    std::string sdbg_file = "toy_example/graph/graph";
    cout << "Loading the graph..." << endl;
    SDBG sdbg;
    sdbg.LoadFromFile(sdbg_file.c_str());
    cout << "Loaded the graph" << endl;
    cout << "Cycle Algorithm Start" << endl;
    CycleFinder cycle_finder(sdbg, 100, 47, "toy_example");
    cout << "Cycle Algorithm End" << endl;
*/
    if (std::string(argv[1])=="u_mode" && argc == 3 ){
          
        //print each line from "welcome.txt"
        std::ifstream welcome_file("/vol/d/development/git/mCAAT/src/welcome.txt");
        std::string line;
        while (std::getline(welcome_file, line))
        {
            std::cout << line << std::endl;
        }
        welcome_file.close();
        //std::cout << "Usage: ./cycle_finder <mode> <genome_name> <length_bound>" << std::endl;
        std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[2])+"/graph/graph";
        int length_bound = 100;
        
        SDBG sdbg;
        cout << "Loading the graph..." << endl;
        //write to file log.txt
        ofstream log_file;
        string log_filename = "log.txt";
        log_file.open(log_filename,std::ios_base::app);
        log_file << "LOADING THE GRAPH..." << endl;
        log_file.close();
        sdbg.LoadFromFile(sdbg_file.c_str());
        log_file.open(log_filename,std::ios_base::app);
        log_file << "LOADED THE GRAPH " <<std::string(argv[2]) << endl;
        log_file << "Cycle Algorithm Start" << endl;
        log_file.close();
        //GraphWriter gw(sdbg);
        //gw.WriteNodes();
        cout << "Loaded the graph\nPress something\n" << endl; 
        std::cin.ignore();

        cout << "Cycle Algorithm Start" << endl;
        //time the cycle finding algorithm
        auto start_time = std::chrono::high_resolution_clock::now();
        CycleFinder cycle_finder(sdbg, length_bound, 47, std::string(argv[2]).c_str());
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        //convert to minutes
        double duration_in_minutes = duration.count() / 60.0;
        //convert to hours
        double duration_in_hours = duration_in_minutes / 60.0;
        //write to file log.txt
        log_file.open(log_filename,std::ios_base::app);
        log_file << "IN HRS: " << duration_in_hours << " HOURS" << endl;
        log_file << "IN MINUTES: " << duration_in_minutes << " MINUTES" << endl;
        log_file << "IN SECONDS: " << duration.count() << " SECONDS" << endl;
        cout << "Cycle Algorithm Finished in " << duration_in_minutes << " minutes" << endl;
        log_file.close();
        cout << "Cycle Algorithm End" << endl;
        
    }
    if (argc ==2)
    {
        //read in from proof_of_concept/data/<genome_name>/cycles/id_paths_no_duplicates.txt
        //get label from sdbg for each node and write to proof_of_concept/data/<genome_name>/cycles/str_true_positives.txt
        std::string genome_name = std::string(argv[1]);
        std::string id_path_file = "proof_of_concept/data/" + genome_name + "/cycles/id_paths.txt";
        std::string str_true_positives_file = "proof_of_concept/data/" + genome_name + "/cycles/string_after_filter.txt";
        std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[1]) + "/graph/graph";
        SDBG sdbg;
        sdbg.LoadFromFile(sdbg_file.c_str());
        std::ifstream id_path_stream(id_path_file);
        std::ofstream str_true_positives_stream(str_true_positives_file);
        std::string line;
        
        while (std::getline(id_path_stream, line))
        {
            std::vector<std::string> path;
            std::stringstream ss(line);
            std::string node_id;
            while (ss >> node_id)
            {
                path.push_back(FetchNodeLabel(sdbg, std::stoull(node_id)));
            }
            for (int i = 0; i < path.size(); i++)
            {
                str_true_positives_stream << path[i] << " ";
            }
            str_true_positives_stream << "\n";
        }
        id_path_stream.close();
        str_true_positives_stream.close();

    }
    if (std::string(argv[1]) == "g_o")
    {
        std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[2])+"/genome/graph/graph";
        SDBG sdbg;
        sdbg.LoadFromFile(sdbg_file.c_str());
        GraphWriter gw(sdbg);
        gw.GenerateGFAFile();
    }
    if (std::string(argv[1]) == "filter")
    {
        std::string results = "proof_of_concept/data/" + std::string(argv[2])+"/cycles_genome/results.txt";
        //FilterOut(results);

        
    }
    /*
    std::string sdbg_file = "proof_of_concept/data/" + std::string(argv[2])+"/genome/graph/graph";
    SDBG sdbg;
    sdbg.LoadFromFile(sdbg_file.c_str());
    GraphAnalysis ga(sdbg);
    ga.Analysis();*/
    return 0;
}
