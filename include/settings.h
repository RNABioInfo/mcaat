//ifndef SETTINGS_H
#ifndef SETTINGS_H

#define SETTINGS_H

#include <iostream>
#include <thread>
#include "hwinfo/hwinfo.h"
#include "hwinfo/utils/unit.h"
#include <sys/stat.h>
#include <map>
#include <thread>
#include <fstream>
#include <iostream>
//std::pair
#include <utility>
using namespace std;

struct Settings {
    
    /// @brief The input files, output folder, RAM in gigabytes and number of threads
    std::string input_files;
    double ram = 0; // In gigabytes
    std::size_t threads= 0;
    std::string graph_folder = "out/graph";
    std::string cycles_folder = "out/cycles";
    std::string output_folder = "out/CRISPR_Arrays.txt";

    //@brief Check number of logical cores
    uint32_t MaxThreads(int num_logical_cores){
        return std::thread::hardware_concurrency();

    }
    

    std::map<std::string, std::pair<bool, std::string>> SettingsMap(double total_in_b, int threads){
        
        std::map<std::string, std::pair<bool, std::string>> settings_message_map;
        //format the ram to 2 decimal places
        std::string ram_to_check = to_string(this->ram).substr(0, to_string(this->ram).find(".") + 3);
        std::string total_ram_string = to_string(total_in_b).substr(0, to_string(total_in_b).find(".") + 3);
        settings_message_map["Input Files"] = make_pair(true, input_files+" exist(s)");
        settings_message_map["RAM"] =this->ram >1 && this->ram <= total_in_b ? make_pair(true, ram_to_check+" GB") : make_pair(false, "value "+ram_to_check+" is more " + "than total "+ total_ram_string+" amount of RAM(in Bytes)");
        settings_message_map["Threads"] = this->threads > 0 && this->threads <= threads ? make_pair(true, to_string(this->threads) + " are used") : make_pair(false, "value "+to_string(this->threads)+" is more than total "+to_string(threads)+" amount of threads" );
        this->ram = this->ram;
        return settings_message_map;
    }
    /// @brief Prints the summary of the settings
    std::string PrintSettings()
    {
        hwinfo::Memory memory;
        double total_in_b = hwinfo::unit::bytes_to_MiB(memory.total_Bytes());
        std::string erroneous_properties = "";
        int num_logical_cores = std::thread::hardware_concurrency();
        std::map<std::string, std::pair<bool, std::string>> settings_message_map = this->SettingsMap(total_in_b, num_logical_cores);
        for (auto it = settings_message_map.begin(); it != settings_message_map.end(); ++it) {
            if (it->second.first) {
                std::cout << "[✔] " << it->first << ": " << it->second.second << std::endl;
            }
            else {
                erroneous_properties += it->first + " ";
                std::cout <<  "[✗] "  << it->first << ": " << it->second.second <<  std::endl;
            }
        }
        return erroneous_properties;
            
    }
};

#endif
