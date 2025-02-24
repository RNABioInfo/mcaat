//test myLib.h
#include <iostream>
#include "sdbg/sdbg.h"
#include "cycle_finder.h"

#include <fstream>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "hwinfo/hwinfo.h"
#include "hwinfo/utils/unit.h"
#include "filters.h"
#include <cstring>
#include "sdbg_build.h"
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
bool check_for_error(Settings& settings){
    cout<<"Step 1. Checking the inputs: "<<endl;
    string erroneous_message = settings.PrintSettings();
    if(erroneous_message.empty()){
        cout << "All inputs are correct. [✔]" << endl;
        return false;
    }
    cout << "Please check the following: " << erroneous_message << endl;
    return true;
    
}

Settings parseArguments(int argc, char* argv[]) {
    vector<string> input_files_default;
    Settings settings;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input_files") {
            while (++i < argc && argv[i][0] != '-') {
                input_files_default.push_back(argv[i]);
            }
            --i;
        } else if (arg == "--ram") {
            if (++i < argc) {
                settings.ram = std::stod(argv[i]);
            } else {
                throw std::runtime_error("Error: Missing value for --ram");
            }
        } else if (arg == "--threads") {
            if (++i < argc) {
                settings.threads = std::stoul(argv[i]);
            } else {
                throw std::runtime_error("Error: Missing value for --threads");
            }
        } else if (arg == "--graph-folder") {
            if (++i < argc) {
                settings.graph_folder = argv[i];
            } else {
                throw std::runtime_error("Error: Missing value for --graph-folder");
            }
        } else if (arg == "--cycles-folder") {
            if (++i < argc) {
                settings.cycles_folder = argv[i];
            } else {
                throw std::runtime_error("Error: Missing value for --cycles-folder");
            }
        } else if (arg == "--output-folder") {
            if (++i < argc) {
                settings.output_folder = argv[i];
            } else {
                throw std::runtime_error("Error: Missing value for --output-folder");
            }
        }
    }

    // Validate input files
    if (input_files_default.size() < 1 || input_files_default.size() > 2) {
        throw std::runtime_error("Error: You must provide one or two input files.");
    }
    for (const auto& file : input_files_default) {
        if (!std::filesystem::exists(file)) {
            throw std::runtime_error("Error: Input file " + file + " does not exist.");
        }
        settings.input_files+=file+" ";
    }
    

    // Set default values if necessary
    if (settings.threads == 0) {
        settings.threads = std::thread::hardware_concurrency()-2;
    }
    
    hwinfo::Memory memory;
    if (settings.ram == 0) 
          settings.ram = hwinfo::unit::bytes_to_MiB(memory.total_Bytes())*0.9;
    return settings;
}

int main(int argc, char** argv) {
    Settings settings = parseArguments(argc, argv);
    string name_of_genome = "test";
    if(check_for_error(settings)) return 1;
    SDBGBuild sdbg_build(settings);
    int length_bound = 77;
    SDBG sdbg;
    // convert it to char * to be able to use it in the function
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph\nPress something\n" << endl;
    
    cout << "Cycle Algorithm Start" << endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    CycleFinder cycle_finder(sdbg, length_bound, 27, std::string(argv[2]).c_str(),settings.threads);
    int number_of_spacers_total = 0;
    auto cycles = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles.size() << endl;
    Filters filters(sdbg,cycles);
    string results_file = settings.output_folder;
    int number_of_spacers = filters.WriteToFile(results_file);
    cout<<"Saved in: "<<results_file<<endl;
    number_of_spacers_total += number_of_spacers;
    cout << "Number of spacers: " << number_of_spacers_total << endl;
            
}
