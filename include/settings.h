#include <iostream>
#include <string>
#include <thread>
#include "hwinfo/hwinfo.h"
#include "hwinfo/utils/unit.h"

struct settings {
    std::string input_files;
    std::string output_folder;
    std::size_t ram; // In gigabytes
    std::size_t threads;

    // Constructor to initialize default values
    settings(string input_files,string output_folder,int ram, int threads) : 
        ram(ram),
        threads(threads),
        output_folder(output_folder),
        input_files(input_files)
    {
        this->physical_ram_gb();
        this->threads = std::thread::hardware_concurrency();
        if (output_folder!="output") this->output_folder = output_folder;

    }

    // Method to get physical RAM in gigabytes (platform-specific)
    void physical_ram_gb() {
        hwinfo::Memory memory;
        double total_in_gb = hwinfo::unit::bytes_to_MiB(memory.total_Bytes()) * 0.001;
        this->ram = total_in_gb;
    }
    
};
