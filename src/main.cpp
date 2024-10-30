//test myLib.h
#include <iostream>
#include "sdbg/sdbg.h"
#include "cycle_finder.h"
#include "system_resources.h"
#include <fstream>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "hwinfo/hwinfo.h"
#include "hwinfo/utils/unit.h"
//#include "filtering_seqan.h"
//include <seqan/align.h>

///using namespace seqan;
using namespace std;
using std::istringstream;
void print_usage(const char* program_name) {
    cout << "-------------------------------------------------------" << endl;
    cout<<"\n";
    cout << "mCAAT - Metagenomic CRISPR Array Analysis Tool v. 0.1" << endl;
    cout<<"\n";
    cout << "-------------------------------------------------------" << endl;
}

string FetchNodeLabel(size_t node, SDBG& sdbg) {
    std::string label;            
    uint8_t seq[sdbg.k()];
    uint32_t t = sdbg.GetLabel(node, seq);
    for (int i = sdbg.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}
int main(int argc, char** argv) {
    print_usage(argv[0]);
    std::string name_of_genome = "mg_04";
    std::string sdbg_file = "/vol/d/development/git/mCAAT/proof_of_concept/data/"+name_of_genome+"/graph/graph";
    if (argc ==3)
    {
        //read in from proof_of_concept/data/<genome_name>/cycles/id_paths_no_duplicates.txt
        //get label from sdbg for each node and write to proof_of_concept/data/<genome_name>/cycles/str_true_positives.txt
        std::string genome_name = std::string(argv[1]);
        std::string id_path_file = "/vol/d/development/git/mCAAT/proof_of_concept/data/" + genome_name + "/cycles/group_cycles.txt";
        std::string str_true_positives_file = "/vol/d/development/git/mCAAT/proof_of_concept/data/" + genome_name + "/cycles/string_after_filter.txt";
        
       
        

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
                path.push_back(FetchNodeLabel(std::stoull(node_id), sdbg));
            }
            for (int i = 0; i < path.size(); i++)
            {
                str_true_positives_stream << path[i] << " ";
            }
            str_true_positives_stream << "\n";
        }
        id_path_stream.close();
        str_true_positives_stream.close();

    }else{
        int length_bound = 100;
        settings settings("1","2",3,4);
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
        //std::cin.ignore();

        cout << "Cycle Algorithm Start" << endl;
        //time the cycle finding algorithm
        auto start_time = std::chrono::high_resolution_clock::now();
        CycleFinder cycle_finder(sdbg, length_bound, 47, name_of_genome);
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
        
}
